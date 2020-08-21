#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

using namespace std;

double RandomDeviceDouble(double min, double max) {
  random_device rd;
  double f = (double)rd() / (double)rd.max();
  return min + f * (max - min);
}

void RandomFinalReport(string file_name, int n, int map_size,
                       int start_snps_from_id = 1,
                       int start_samples_from_id = 1) {
  random_device rd;
  char alleles[5] = "ATCG", ab[3] = "AB";
  ofstream outfile;
  outfile.open(file_name, ofstream::out | ofstream::trunc);

  outfile << "[Header]\n[Data]\n";
  outfile
      << "SNP Name\tSample ID\tAllele1 - Forward\tAllele2 - Forward\tAllele1 - "
         "Top\tAllele2 - Top\tAllele1 - AB\tAllele2 - AB\tGC Score\tX\tY\n";
  for (int i = 0; i < n * map_size; i++) {
    outfile << "SNP" << to_string(start_snps_from_id + i % map_size) << "\t"
            << "SAM" << to_string(start_samples_from_id + i / map_size) << "\t"
            << alleles[rd() % 4] << "\t" << alleles[rd() % 4] << "\t"
            << alleles[rd() % 4] << "\t" << alleles[rd() % 4] << "\t"
            << ab[rd() % 2] << "\t" << ab[rd() % 2] << "\t" << setprecision(4)
            << RandomDeviceDouble(0.0, 1.0) << "\t" << setprecision(3)
            << RandomDeviceDouble(0.0, 1.0) << "\t" << setprecision(3)
            << RandomDeviceDouble(0.0, 1.0) << "\n";
  }

  outfile.close();
}

// Main function
// Arguments:
//  argc ==
//  argv = [file_name, number_of_samples, map_size,
//          start_snps_from_id, start_samples_from_id]
int main(int argc, char const *argv[]) {
  // get parameters from command-line
  if (argc < 4 || argc > 6)
    throw std::invalid_argument(
        "usage: main file_name number_of_samples map_size [start_snps_from_id "
        "start_samples_from_id]");

  string file_name = argv[1];
  int n = atoi(argv[2]), map_size = atoi(argv[3]);
  // int seed = time(nullptr);
  int start_snps_from_id = 1, start_samples_from_id = 1;
  // if (argc > 4) seed = atoi(argv[4]);
  if (argc > 4) start_snps_from_id = atoi(argv[4]);
  if (argc > 5) start_samples_from_id = atoi(argv[5]);

  cout << "Generating Final Report File (" << file_name << "): " << n
       << " samples, " << map_size << "SNPs" << endl;

  RandomFinalReport(file_name, n, map_size, start_snps_from_id,
                    start_samples_from_id);
  auto file_size =
      ifstream(file_name, ifstream::ate | ifstream::binary).tellg();
  cout << "Generated Final Report File: " << (double)file_size / (1024 * 1024)
       << "MB" << endl;

  return 0;
}
