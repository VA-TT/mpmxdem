#include "Core/MPMbox.hpp"
#include "Core/MaterialPoint.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <string>
#include <vector>

#include "ConstitutiveModels/CHCL_DEM.hpp"
#include "ConstitutiveModels/CHCL_DEM_Multi.hpp"
#include "ConstitutiveModels/MohrCoulomb.hpp"

namespace fs = std::filesystem;

struct Options {
  fs::path directory  = ".";
  std::string mode    = "smooth";
  int bins            = 80;
  int minimumBinCount = 2;
  float fitIMin       = 1.0e-4f;
  float fitIMax       = 1.0e-1f;
};

struct Row {
  float inertialNumber;
  float mohrFriction;
  float rheologyFriction;
  float packingFraction;
};

struct Bin {
  int count              = 0;
  double sumI            = 0.0;
  double sumValue        = 0.0;
  double sumValueSquared = 0.0;
};

struct RheologyFile {
  int number;
  fs::path path;
};

void applyRawMaterialPointFields(const MPMbox &box, std::vector<ProcessedDataMP> &data) {
  for (size_t index = 0; index < box.MP.size(); ++index) {
    data[index].stress           = box.MP[index].stress;
    data[index].outOfPlaneStress = box.MP[index].outOfPlaneStress;
    data[index].rho              = box.MP[index].density;
    data[index].outOfPlaneShearXZ = box.MP[index].outOfPlaneShearXZ;
    data[index].outOfPlaneShearYZ = box.MP[index].outOfPlaneShearYZ;
  }
}

void printUsage(const char *program) {
  std::cerr << "Usage: " << program
            << " DIRECTORY [--mode smooth|raw] [--bins N] [--min-count N] [--fit-min I] [--fit-max I]\n";
}

Options parseOptions(int argc, char **argv) {
  if (argc < 2) {
    printUsage(argv[0]);
    throw std::runtime_error("missing input directory");
  }

  Options options;
  options.directory = argv[1];
  for (int i = 2; i < argc; ++i) {
    const std::string argument = argv[i];
    auto readValue             = [&](const char *name) {
      if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
      return std::string(argv[++i]);
    };
    if (argument == "--mode") {
      options.mode = readValue("--mode");
      if (options.mode != "smooth" && options.mode != "raw") {
        throw std::runtime_error("--mode must be smooth or raw");
      }
    } else if (argument == "--bins") {
      options.bins = std::stoi(readValue("--bins"));
    } else if (argument == "--min-count") {
      options.minimumBinCount = std::stoi(readValue("--min-count"));
    } else if (argument == "--fit-min") {
      options.fitIMin = std::stof(readValue("--fit-min"));
    } else if (argument == "--fit-max") {
      options.fitIMax = std::stof(readValue("--fit-max"));
    } else {
      throw std::runtime_error("unknown option: " + argument);
    }
  }
  if (options.bins <= 0 || options.minimumBinCount <= 0 || options.fitIMin <= 0.0f ||
      options.fitIMax <= options.fitIMin) {
    throw std::runtime_error("invalid binning or fit options");
  }
  return options;
}

std::vector<RheologyFile> findConfigurationFiles(const fs::path &directory) {
  const std::regex pattern(R"(^conf([0-9]+)\.txt$)");
  std::vector<RheologyFile> files;
  for (const auto &entry : fs::directory_iterator(directory)) {
    if (!entry.is_regular_file()) continue;
    std::smatch match;
    const std::string name = entry.path().filename().string();
    if (std::regex_match(name, match, pattern)) files.push_back({std::stoi(match[1]), entry.path()});
  }
  std::sort(files.begin(), files.end(),
            [](const RheologyFile &left, const RheologyFile &right) { return left.number < right.number; });
  return files;
}

float computeShearRate(const ProcessedDataMP &data) {
  const float normalDifference = static_cast<float>(data.velGrad.xx - data.velGrad.yy);
  const float shearSum         = static_cast<float>(data.velGrad.xy + data.velGrad.yx);
  return std::sqrt(normalDifference * normalDifference + shearSum * shearSum);
}

float computePressure(const MPMbox &box, const std::vector<ProcessedDataMP> &data, size_t index) {
  const float sx = static_cast<float>(data[index].stress.xx);
  const float sy = static_cast<float>(data[index].stress.yy);
  const float sz = box.MP[index].isDoubleScale ? static_cast<float>(data[index].outOfPlaneStress) : 0.0f;
  return box.MP[index].isDoubleScale ? -(sx + sy + sz) / 3.0f : -(sx + sy) / 2.0f;
}

float representativeDiameter(const MPMbox &box, size_t index, std::map<std::string, float> &diameters) {
  const auto *multi = dynamic_cast<const CHCL_DEM_Multi *>(box.MP[index].constitutiveModel);
  const auto *chcl  = dynamic_cast<const CHCL_DEM *>(box.MP[index].constitutiveModel);
  if (multi == nullptr && chcl == nullptr) return 0.0f;

  const std::string fileName = multi != nullptr ? multi->getFileNameFor(box.MP[index]) : chcl->fileName;
  const auto cached          = diameters.find(fileName);
  if (cached != diameters.end()) return cached->second;

  PBC3Dbox sample;
  sample.loadConf(fileName.c_str());
  const float diameter = 2.0f * static_cast<float>(sample.Rmean);
  diameters[fileName]  = diameter;
  return diameter;
}

float grainDensity(const MPMbox &box, size_t index, std::map<std::string, float> &densities) {
  const auto *multi = dynamic_cast<const CHCL_DEM_Multi *>(box.MP[index].constitutiveModel);
  const auto *chcl  = dynamic_cast<const CHCL_DEM *>(box.MP[index].constitutiveModel);
  if (multi == nullptr && chcl == nullptr) return 0.0f;

  const std::string fileName = multi != nullptr ? multi->getFileNameFor(box.MP[index]) : chcl->fileName;
  const auto cached          = densities.find(fileName);
  if (cached != densities.end()) return cached->second;

  PBC3Dbox sample;
  sample.loadConf(fileName.c_str());
  const float density = static_cast<float>(sample.density);
  densities[fileName] = density;
  return density;
}

float computeMohrFriction(const MPMbox &box, const std::vector<ProcessedDataMP> &data, size_t index) {
  const float sx         = static_cast<float>(data[index].stress.xx);
  const float sy         = static_cast<float>(data[index].stress.yy);
  const float sxy        = 0.5f * static_cast<float>(data[index].stress.xy + data[index].stress.yx);
  const float sz         = box.MP[index].isDoubleScale ? static_cast<float>(data[index].outOfPlaneStress) : 0.0f;
  const float radius2D   = std::sqrt(0.25f * (sx - sy) * (sx - sy) + sxy * sxy);
  const float sigmaA     = -0.5f * (sx + sy) + radius2D;
  const float sigmaB     = -0.5f * (sx + sy) - radius2D;
  const float sigma1     = box.MP[index].isDoubleScale ? std::max({sigmaA, sigmaB, -sz}) : std::max(sigmaA, sigmaB);
  const float sigma3     = box.MP[index].isDoubleScale ? std::min({sigmaA, sigmaB, -sz}) : std::min(sigmaA, sigmaB);
  const float mohrCenter = 0.5f * (sigma1 + sigma3);
  const float mohrRadius = 0.5f * (sigma1 - sigma3);
  const auto *model      = dynamic_cast<const MohrCoulomb *>(box.MP[index].constitutiveModel);
  const float cohesion   = model == nullptr ? 0.0f : std::max(0.0f, static_cast<float>(model->Cohesion));
  const float amplitude  = std::hypot(mohrCenter, cohesion);
  if (amplitude <= 1.0e-10f || mohrRadius >= amplitude) return 0.0f;
  const float phi = std::asin(mohrRadius / amplitude) - std::atan2(cohesion, mohrCenter);
  return std::sin(std::max(0.0f, phi));
}

Row computeRow(const MPMbox &box, const std::vector<ProcessedDataMP> &data, size_t index,
               std::map<std::string, float> &diameters, std::map<std::string, float> &densities) {
  const float pressure  = computePressure(box, data, index);
  const float shearRate = computeShearRate(data[index]);
  const float mass      = static_cast<float>(data[index].rho * box.MP[index].vol);
  float inertialNumber  = 0.0f;
  if (box.MP[index].isDoubleScale) {
    const float diameter = representativeDiameter(box, index, diameters);
    if (diameter > 0.0f) inertialNumber = shearRate * std::sqrt(mass / (diameter * (std::fabs(pressure) + 1.0e-6f)));
  } else {
    inertialNumber = shearRate * std::sqrt(mass / (std::fabs(pressure) + 1.0e-6f));
  }

  const float sx               = static_cast<float>(data[index].stress.xx);
  const float sy               = static_cast<float>(data[index].stress.yy);
  const float sxy              = 0.5f * static_cast<float>(data[index].stress.xy + data[index].stress.yx);
  const float sz               = box.MP[index].isDoubleScale ? static_cast<float>(data[index].outOfPlaneStress) : 0.0f;
  const float sxx              = sx + pressure;
  const float syy              = sy + pressure;
  const float szz              = box.MP[index].isDoubleScale ? sz + pressure : 0.0f;
  const float sxz = box.MP[index].isDoubleScale ? static_cast<float>(data[index].outOfPlaneShearXZ) : 0.0f;
  const float syz = box.MP[index].isDoubleScale ? static_cast<float>(data[index].outOfPlaneShearYZ) : 0.0f;
  const float tau = std::sqrt(0.5f * (sxx * sxx + syy * syy + szz * szz + 2.0f * (sxy * sxy + sxz * sxz + syz * syz)));
  const float rheologyFriction = pressure > 1.0e-6f ? tau / pressure : 0.0f;
  const float density          = grainDensity(box, index, densities);
  const float packingFraction  = density > 0.0f ? static_cast<float>(data[index].rho) / density : 0.0f;
  return {inertialNumber, computeMohrFriction(box, data, index), rheologyFriction, packingFraction};
}

void writeRheologyFile(const fs::path &path, const MPMbox &box, const std::vector<ProcessedDataMP> &data,
                       std::vector<Row> &rows, const std::string &mode) {
  std::ofstream output(path);
  output << "# MP_index MP_x MP_y I mu_mohr_sin_phi mu_rheology_tau_over_p packing_fraction_Phi (" << mode
         << " stress/density, smoothed velocity gradient)\n";
  std::map<std::string, float> diameters;
  std::map<std::string, float> densities;
  rows.resize(box.MP.size());
  output << std::setprecision(12);
  for (size_t index = 0; index < box.MP.size(); ++index) {
    rows[index] = computeRow(box, data, index, diameters, densities);
    output << index << ' ' << box.MP[index].pos.x << ' ' << box.MP[index].pos.y << ' ' << rows[index].inertialNumber
           << ' ' << rows[index].mohrFriction << ' ' << rows[index].rheologyFriction << ' '
           << rows[index].packingFraction << '\n';
  }
}

void addToBin(Bin &bin, float inertialNumber, float value) {
  ++bin.count;
  bin.sumI += inertialNumber;
  bin.sumValue += value;
  bin.sumValueSquared += static_cast<double>(value) * value;
}

void writeBinnedFile(const fs::path &path, const std::vector<Row> &rows, int bins, int minimumCount, float valueMin,
                     float valueMax, float binIMin, float binIMax, const char *valueName) {
  std::vector<Bin> values(static_cast<size_t>(bins));
  const double logMin   = std::log10(binIMin);
  const double logWidth = (std::log10(binIMax) - logMin) / bins;
  for (const Row &row : rows) {
    const float value = std::string(valueName) == "mu" ? row.rheologyFriction : row.packingFraction;
    if (row.inertialNumber < binIMin || row.inertialNumber > binIMax || value < valueMin || value > valueMax) {
      continue;
    }
    int bin = static_cast<int>((std::log10(row.inertialNumber) - logMin) / logWidth);
    bin     = std::clamp(bin, 0, bins - 1);
    addToBin(values[static_cast<size_t>(bin)], row.inertialNumber, value);
  }

  std::ofstream output(path);
  output << "# I_mean " << (std::string(valueName) == "mu" ? "mu_mean" : "Phi_mean") << " "
         << (std::string(valueName) == "mu" ? "mu_std" : "Phi_std") << " n I_left I_right\n";
  output << std::setprecision(12);
  for (int index = 0; index < bins; ++index) {
    const Bin &bin = values[static_cast<size_t>(index)];
    if (bin.count < minimumCount) continue;
    const double mean     = bin.sumValue / bin.count;
    const double variance = std::max(0.0, bin.sumValueSquared / bin.count - mean * mean);
    const double left     = std::pow(10.0, logMin + index * logWidth);
    const double right    = std::pow(10.0, logMin + (index + 1) * logWidth);
    output << bin.sumI / bin.count << ' ' << mean << ' ' << std::sqrt(variance) << ' ' << bin.count << ' ' << left
           << ' ' << right << '\n';
  }
}

int main(int argc, char **argv) {
  INIT_TIMERS();
  try {
    const Options options     = parseOptions(argc, argv);
    const fs::path directory  = fs::absolute(options.directory);
    const auto configurations = findConfigurationFiles(directory);
    if (configurations.empty()) throw std::runtime_error("no confN.txt files found");

    const fs::path originalDirectory = fs::current_path();
    fs::current_path(directory);
    std::ofstream meanOutput(directory / "mean_I_time.dat");
    meanOutput << "# conf t I_mean I_std n file\n" << std::setprecision(12);
    std::vector<Row> allRows;
    float globalIMin = std::numeric_limits<float>::max();
    float globalIMax = 0.0f;

    for (const RheologyFile &configuration : configurations) {
      MPMbox box;
      box.computationMode = false;
      box.read(configuration.path.filename().string().c_str());
      std::vector<ProcessedDataMP> smoothedData;
      box.postProcess(smoothedData);
      if (options.mode == "raw") applyRawMaterialPointFields(box, smoothedData);
      std::vector<Row> rows;
      const fs::path outputName = directory / ("rheology_conf" + std::to_string(configuration.number) + ".txt");
      writeRheologyFile(outputName, box, smoothedData, rows, options.mode);
      allRows.insert(allRows.end(), rows.begin(), rows.end());

      double sumI        = 0.0;
      double sumISquared = 0.0;
      int count          = 0;
      for (const Row &row : rows) {
        if (row.inertialNumber <= 0.0f) continue;
        sumI += row.inertialNumber;
        sumISquared += static_cast<double>(row.inertialNumber) * row.inertialNumber;
        globalIMin = std::min(globalIMin, row.inertialNumber);
        globalIMax = std::max(globalIMax, row.inertialNumber);
        ++count;
      }
      if (count > 0) {
        const double mean     = sumI / count;
        const double variance = std::max(0.0, sumISquared / count - mean * mean);
        meanOutput << configuration.number << ' ' << box.t << ' ' << mean << ' ' << std::sqrt(variance) << ' ' << count
                   << ' ' << outputName.filename().string() << '\n';
      }
      std::cout << "Created " << outputName << '\n';
    }
    fs::current_path(originalDirectory);

    if (globalIMin < globalIMax) {
      const float binIMin = std::max(options.fitIMin, globalIMin);
      const float binIMax = std::min(options.fitIMax * 10.0f, globalIMax);
      writeBinnedFile(directory / "binned_mu_stats.dat", allRows, options.bins, options.minimumBinCount, 0.0f, 1.0f,
                      binIMin, binIMax, "mu");
      writeBinnedFile(directory / "binned_packing_stats.dat", allRows, options.bins, options.minimumBinCount, 0.0f,
                      1.0f, binIMin, binIMax, "Phi");
    }
    std::cout << "Created " << directory / "mean_I_time.dat" << '\n';
  } catch (const std::exception &error) {
    std::cerr << "rheology: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
