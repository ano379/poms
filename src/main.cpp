#include "test.hpp"
#include "syn_data.hpp"
using namespace std;
int main()
{
    bool test_poms_prob = true;
    bool test_hpdfs_time = true;
    bool test_poms_time = true;
    bool get_level_stats = true;

// ==========================================================================
    if(test_poms_prob){
        // poms prob
        string f_name = "../data/amazon.txt";
        string of_name = "../result/amazon_poms_k";
        vector<int> ks = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        test_poms_k(f_name, f_name+"_query", of_name, of_name+"_stats", ks);
        f_name = "../data/wordnet.txt";
        of_name = "../result/wordnet_poms_k";
        test_poms_k(f_name, f_name+"_query", of_name, of_name+"_stats", ks);
        f_name = "../data/imagenet.txt";
        of_name = "../result/imagenet_poms_k";
        test_poms_k(f_name, f_name+"_query", of_name, of_name+"_stats", ks);
        int n = 1000000, d = 30;
        float r = 0.1;
        int k = 5;
        int sample = 1000;
        poms_vary_n(k, d, r);
        poms_vary_d(n, k, r);
        poms_vary_r(n, k, d);
        f_name = "../data/syn/syn_" + to_string(n/10000) + "w_d" + to_string(d) + "_r0" + to_string((int)(r*10));
        of_name = "../result/syn_poms_vary_k";
        poms_vary_k(f_name, f_name+"_query", of_name, of_name+"_stats");
    }
// ==========================================================================
if(test_hpdfs_time){
        string f_name = "../data/amazon.txt";
        string of_name = "../result/real_hpdfs_time";
        ofstream ofile(of_name);
        ofile<<"file, naive_time, nm_time, bridge_time, nm_speedup, bridge_speedup"<<endl;
        ofile.close();
        test_hpdfs_f(f_name, of_name);
        f_name = "../data/wordnet.txt";
        test_hpdfs_f(f_name, of_name);
        f_name = "../data/imagenet.txt";
        test_hpdfs_f(f_name, of_name);
        test_hpdfs_syn();
    }
// ==========================================================================
if(test_poms_time){
        string f_name = "../data/amazon.txt";
        string of_name = "../result/real_poms_time";
        int sample_size = 1000; 
        // int sample_size = 0; //for testing all the leaves
        poms_time(f_name, of_name, sample_size);
        f_name = "../data/wordnet.txt";
        poms_time(f_name, of_name, sample_size);
        f_name = "../data/imagenet.txt";
        poms_time(f_name, of_name, sample_size);
        int n = 1000000, d = 30;
        float r = 0.1;
        int sample = 1000;
        of_name = "../result/syn_poms_vary_n_time";
        poms_vary_n_time(of_name, d, r, sample);
        of_name = "../result/syn_poms_vary_d_time";
        poms_vary_d_time(of_name, n , r, sample);
        of_name = "../result/syn_poms_vary_r_time";
        poms_vary_r_time(of_name, n, d, sample);
}
// =============================================================================
if(get_level_stats){
        string f_name = "../data/amazon.txt";
        string of_name = "../result/amazon_poms_k";
        vector<int> ks = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        string level_file = "../result/amazon_poms_lvl";
        process_lvl(f_name, of_name+"_stats", f_name+"_query", level_file, ks);
        f_name = "../data/wordnet.txt";
        of_name = "../result/wordnet_poms_k";
        level_file = "../result/wordnet_poms_lvl";
        process_lvl(f_name, of_name+"_stats", f_name+"_query", level_file, ks);
        f_name = "../data/imagenet.txt";
        of_name = "../result/imagenet_poms_k";
        level_file = "../result/imagenet_poms_lvl";
        process_lvl(f_name, of_name+"_stats", f_name+"_query", level_file, ks);
}
    return 0;
}
