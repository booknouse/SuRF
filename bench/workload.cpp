#include "bench.hpp"
#include "filter_factory.hpp"

int main(int argc, char *argv[]) {
    if (argc != 9) {
	std::cout << "Usage:\n";
	std::cout << "1. filter type: SuRF, SuRFHash, SuRFReal, Bloom, ARF\n";
	std::cout << "2. workload type: mixed, alterByte (only for email key)\n";
	std::cout << "3. percentage of keys inserted: 0 < num <= 100\n";
	std::cout << "4. byte position (conting from last, only for alterByte): num\n";
	std::cout << "5. key type: randint, timestamp, email\n";
	std::cout << "6. query type: point, range\n";
	std::cout << "7. range size: num\n";
	std::cout << "8. distribution: uniform, zipfian, latest\n";
	return -1;
    }

    std::string filter_type = argv[1];
    std::string workload_type = argv[2];
    unsigned percent = atoi(argv[3]);
    unsigned byte_pos = atoi(argv[4]);
    std::string key_type = argv[5];
    std::string query_type = argv[6];
    uint64_t range_size = atoi(argv[7]);
    std::string distribution = argv[8];

    // check args ====================================================
    if (filter_type.compare(std::string("SuRF")) != 0
	&& filter_type.compare(std::string("SuRFHash")) != 0
	&& filter_type.compare(std::string("SuRFReal")) != 0
	&& filter_type.compare(std::string("Bloom")) != 0
	&& filter_type.compare(std::string("ARF")) != 0) {
	std::cout << bench::kRed << "WRONG filter type\n" << bench::kNoColor;
	return -1;
    }

    if (workload_type.compare(std::string("mixed")) != 0
	&& workload_type.compare(std::string("alterByte")) == 0) {
	std::cout << bench::kRed << "WRONG workload type\n" << bench::kNoColor;
	return -1;
    }

    if (percent > 100) {
	std::cout << bench::kRed << "WRONG percentage\n" << bench::kNoColor;
	return -1;
    }

    if (key_type.compare(std::string("randint")) != 0
	&& key_type.compare(std::string("timestamp")) != 0
	&& key_type.compare(std::string("email")) != 0) {
	std::cout << bench::kRed << "WRONG key type\n" << bench::kNoColor;
	return -1;
    }

    if (query_type.compare(std::string("point")) != 0
	&& query_type.compare(std::string("range")) != 0) {
	std::cout << bench::kRed << "WRONG query type\n" << bench::kNoColor;
	return -1;
    }

    if (distribution.compare(std::string("uniform")) != 0
	&& distribution.compare(std::string("zipfian")) != 0
	&& distribution.compare(std::string("latest")) != 0) {
	std::cout << bench::kRed << "WRONG distribution\n" << bench::kNoColor;
	return -1;
    }

    // load keys from files =======================================
    std::string load_file = "workloads/load_";
    load_file += key_type;
    std::vector<std::string> load_keys;
    if (key_type.compare(std::string("email")) == 0)
	bench::loadKeysFromFile(load_file, false, load_keys);
    else
	bench::loadKeysFromFile(load_file, true, load_keys);

    std::string txn_file = "workloads/txn_";
    txn_file += key_type;
    txn_file += "_";
    txn_file += distribution;
    std::vector<std::string> txn_keys;
    if (key_type.compare(std::string("email")) == 0)
	bench::loadKeysFromFile(txn_file, false, txn_keys);
    else
	bench::loadKeysFromFile(txn_file, true, txn_keys);

    std::vector<std::string> insert_keys;
    bench::selectKeysToInsert(percent, insert_keys, load_keys);

    if (workload_type.compare(std::string("alterByte")) == 0)
	bench::modifyKeyByte(txn_keys, byte_pos);

    // compute upperbound keys for range queries =================
    std::vector<std::string> upper_bound_keys;
    if (query_type.compare(std::string("range")) == 0) {
	for (int i = 0; i < txn_keys.size(); i++)
	    upper_bound_keys.push_back(bench::getUpperBoundKey(key_type, txn_keys[i], range_size));
    }

    // create filter ==============================================
    bench::Filter* filter = bench::FilterFactory::createFilter(filter_type, insert_keys);

    // execute transactions =======================================
    int64_t positives = 0;
    double start_time = bench::getNow();
    if (query_type.compare(std::string("point")) == 0) {
	for (int i = 0; i < txn_keys.size(); i++)
	    positives += (int)filter->lookup(txn_keys[i]);
    } else if (query_type.compare(std::string("range")) == 0) {
	for (int i = 0; i < txn_keys.size(); i++)
	//for (int i = 0; i < 5; i++)
	    positives += (int)filter->lookupRange(txn_keys[i], upper_bound_keys[i]);
    }
    double end_time = bench::getNow();

    // compute true positives ======================================
    std::map<std::string, bool> ht;
    for (int i = 0; i < insert_keys.size(); i++)
	ht[insert_keys[i]] = true;

    int64_t true_positives = 0;
    std::map<std::string, bool>::iterator ht_iter;
    if (query_type.compare(std::string("point")) == 0) {
	for (int i = 0; i < txn_keys.size(); i++) {
	    ht_iter = ht.find(txn_keys[i]);
	    true_positives += (ht_iter != ht.end());
	}
    } else if (query_type.compare(std::string("range")) == 0) {
	for (int i = 0; i < txn_keys.size(); i++) {
	    ht_iter = ht.upper_bound(txn_keys[i]);
	    if (ht_iter != ht.end()) {
		std::string fetched_key = ht_iter->first;
		true_positives += (fetched_key.compare(upper_bound_keys[i]) < 0);
	    }
	}
    }
    int64_t false_positives = positives - true_positives;
    assert(false_positives >= 0);
    int64_t true_negatives = txn_keys.size() - true_positives;

    // print
    double tput = txn_keys.size() / (end_time - start_time) / 1000000; // Mops/sec
    std::cout << bench::kGreen << "Throughput = " << bench::kNoColor << tput << "\n";

    std::cout << "positives = " << positives << "\n";
    std::cout << "true positives = " << true_positives << "\n";
    std::cout << "false positives = " << false_positives << "\n";
    std::cout << "true negatives = " << true_negatives << "\n";

    double fp_rate = 0;
    if (false_positives > 0)
	fp_rate = false_positives / (true_negatives + false_positives + 0.0);
    std::cout << bench::kGreen << "False Positive Rate = " << bench::kNoColor << fp_rate << "\n";

    std::cout << bench::kGreen << "Memory = " << bench::kNoColor << filter->getMemoryUsage() << "\n\n";

    return 0;
}
