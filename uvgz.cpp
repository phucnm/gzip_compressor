#include <iostream>
#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include "output_stream.hpp"
#include <deque>
#include <algorithm>
#include <queue>
#include <numeric>

#define CRCPP_USE_CPP11
#include "CRC.h"

#define MAX_MATCH_LENGTH 258
#define CHUNK_SIZE 65535
#define WINDOW_SIZE 32768
// Maximum number of hash-chain hit, to avoid the worst case where the chain is entire input
// e.g. "aaa....aaa"
#define MAX_HASH_QUERY 512

using namespace std;

vector<u8> history;

array<pair<u32, u32>, 288> fixed_ll_codes;
array<array<u32, 4>, 29> length_code_ranges {{
    {257,0,3,3},     {258,0,4,4},     {259,0,5,5},     {260,0,6,6},     {261,0,7,7},
    {262,0,8,8},     {263,0,9,9},     {264,0,10,10},   {265,1,11,12},   {266,1,13,14},
    {267,1,15,16},   {268,1,17,18},   {269,2,19,22},   {270,2,23,26},   {271,2,27,30},
    {272,2,31,34},   {273,3,35,42},   {274,3,43,50},   {275,3,51,58},   {276,3,59,66},
    {277,4,67,82},   {278,4,83,98},   {279,4,99,114},  {280,4,115,130}, {281,5,131,162},
    {282,5,163,194}, {283,5,195,226}, {284,5,227,257}, {285,0,258,258}
}};
array<array<u32, 4>, 30> dist_code_ranges {{
    {0,0,1,1},         {1,0,2,2},          {2,0,3,3},           {3,0,4,4},           {4,1,5,6},
    {5,1,7,8},         {6,2,9,12},         {7,2,13,16},         {8,3,17,24},         {9,3,25,32},
    {10,4,33,48},      {11,4,49,64},       {12,5,65,96},        {13,5,97,128},       {14,6,129,192},
    {15,6,193,256},    {16,7,257,384},     {17,7,385,512},      {18,8,513,768},      {19,8,769,1024},
    {20,9,1025,1536},  {21,9,1537,2048},   {22,10,2049,3072},   {23,10,3073,4096},   {24,11,4097,6144},
    {25,11,6145,8192}, {26,12,8193,12288}, {27,12,12289,16384}, {28,13,16385,24576}, {29,13,24577,32768},
}};

array<u16, 19> cl_order {
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
};


/*
 Helper method to reverse bit in the same num_bits range
 */
u32 reverse_bits(u32 val, int num_bits) {
    u32 res = 0;
    u32 _val = val;
    int remain = num_bits;
    while (_val) {
        res <<= 1;
        res |= _val & 1;
        _val >>= 1;
        remain--;
    }
    res <<= remain;
    return res;
}

/*
 Hash-Chain Match Finder
 */

uint32_t hash32(uint32_t x) {
  x ^= x >> 18;
  x *= uint32_t(0xa136aaad);
  x ^= x >> 16;
  x *= uint32_t(0x9f6d62d7);
  x ^= x >> 17;
  return x;
}

class HashChainMatcher {
public:
    HashChainMatcher() {
        this->hash_table.resize(WINDOW_SIZE, -WINDOW_SIZE - 1);
        this->chain.resize(WINDOW_SIZE, -WINDOW_SIZE - 1);
    }
    void reset() {
        this->hash_table.resize(WINDOW_SIZE, -WINDOW_SIZE - 1);
        this->chain.resize(WINDOW_SIZE, -WINDOW_SIZE - 1);
    }
    void insert(const array<u8, CHUNK_SIZE>& bytes, u32 pos) {
        int key = hash32((bytes[pos]) | (bytes[pos+1] << 8) | (bytes[pos+2] << 16)) & (WINDOW_SIZE-1);
        this->chain[pos & (WINDOW_SIZE-1)] = hash_table[key];
        this->hash_table[key] = pos;
    }
    pair<u32, u32> find(const array<u8, CHUNK_SIZE>& block, u32 block_size, int pos) {
        if (pos > block_size - 3) {
            return make_pair(0, 0);
        }
        int key = hash32((block[pos]) | (block[pos+1]<<8) | (block[pos+2]<<16)) & (WINDOW_SIZE-1);

        int match_pos = hash_table[key];
        int count = 1;
        int max_length = 0;
        int dist = 0;
        while (match_pos > pos - WINDOW_SIZE && count++ < MAX_HASH_QUERY) {
            if (match_pos < pos) {
                int len = longest_match(block, block_size, match_pos, pos);
                if (len > max_length) {
                    max_length = len;
                    dist = pos - match_pos;
                }
            }
            match_pos = chain[match_pos & (WINDOW_SIZE-1)];
        }

        chain[pos & (WINDOW_SIZE-1) ] = hash_table[key];
        hash_table[key] = pos;

        return make_pair(max_length, dist);
    }
private:
    int longest_match(const array<u8, CHUNK_SIZE>& block, u32 block_size, int match_pos, int cur_pos) {
        int len = 0;
        while (cur_pos < block_size && block[match_pos] == block[cur_pos] && len < MAX_MATCH_LENGTH) {
            match_pos++;
            cur_pos++;
            len++;
        }
        return len;
    }
    vector<int> hash_table;
    vector<int> chain;
};

HashChainMatcher matcher;

/*
 Huffman Tree Node
 */
class Node {
public:
    Node(u32 symbol, u32 value) {
        this->symbol = symbol;
        this->value = value;
    }

    u32 symbol;
    u32 value;
    Node* left;
    Node* right;
};

/*
 Huffman Tree builder.
 Provides a canonical codes builder and package-merge builder as static methods.
 */

class Compare
{
public:
     bool operator ()(const Node* lhs, const Node* rhs) const
       {
           return lhs->value > rhs->value;
       }
};

class HuffmanCodingBuilder {
public:
    static unordered_map<u32, u32> package_merge(const vector<u32>& freqs, int max_bits) {
        vector<pair<vector<u32>, u32>> original_freqs;
        for (u32 i = 0; i < freqs.size(); i++) {
            if (freqs[i] > 0) {
                original_freqs.push_back(make_pair(vector<u32>{ i}, freqs[i]));
            }
        }
        if (original_freqs.empty()) {
            return unordered_map<u32, u32>{};
        } else if (original_freqs.size() == 1) {
            auto single_sym = original_freqs[0];
            u32 sym = single_sym.first[0];
            return unordered_map<u32, u32>({ {sym, 1} });
        }
        vector<pair<u32, u32>> res;
        vector<pair<vector<u32>, u32>> pkgs(original_freqs);
        sort(pkgs.begin(), pkgs.end(), [](const pair<vector<u32>, u32>& a, const pair<vector<u32>, u32>& b) -> bool {
            return a.second < b.second;
        });
        u32 num_pick_items = (u32)(2 * original_freqs.size() - 2);
        while (--max_bits) {
            if (pkgs.size() % 2 == 1) {
                pkgs.pop_back();
            }
            vector<pair<vector<u32>, u32>> new_pkgs;
            for (size_t i = 0; i < pkgs.size(); i+= 2) {
                vector<u32> merged_sym(pkgs[i].first);
                merged_sym.insert(merged_sym.end(), pkgs[i+1].first.begin(), pkgs[i+1].first.end());

                u32 merge_val = pkgs[i].second + pkgs[i+1].second;
                new_pkgs.push_back(make_pair(merged_sym, merge_val));
            }
            pkgs = vector<pair<vector<u32>, u32>>(original_freqs);
            pkgs.insert(pkgs.begin(), new_pkgs.begin(), new_pkgs.end());
            sort(pkgs.begin(), pkgs.end(), [](const pair<vector<u32>, u32>& a, const pair<vector<u32>, u32>& b) -> bool {
                return a.second < b.second;
            });
        }

        unordered_map<u32, u32> map;
        for (size_t i = 0; i < num_pick_items; i++) {
            auto v = pkgs[i].first;
            for (u32 sym: v) {
                if (map.find(sym) == map.end()) {
                    map[sym] = 1;
                } else {
                    map[sym]++;
                }
            }
        }

        return map;
    }
    static vector<u32> generate_from_freqs(vector<u32> freqs) {
        vector<u32> res(freqs.size());

        //build priority queue
        priority_queue<Node*, vector<Node*>, Compare> q;
        for (u32 i = 0; i < freqs.size(); i++) {
            if (freqs[i] > 0) {
                Node* node = new Node(i, freqs[i]);
                node->left = nullptr;
                node->right = nullptr;
                q.push(node);
            }
        }

        Node *root = nullptr;

        if (q.size() == 1) {
            auto node = q.top();
            res[node->symbol] = 1;
            return res;
        }

        while (q.size() > 1) {
            Node *x = q.top();
            q.pop();
            Node * y = q.top();
            q.pop();

            Node* newNode = new Node(-1, x->value + y->value);
            newNode->left = x;
            newNode->right = y;
            root = newNode;
            q.push(newNode);
        }

        tree_traverse(root, res);

        return res;
    }

    static void tree_traverse(Node *root, vector<u32>& result) {
        do_tree_traverse(root, result, 0);
    }

    static void do_tree_traverse(Node *root, vector<u32>& result, u32 code_length) {
        if (root == nullptr) {
            return;
        }
        if (root->left == nullptr && root->right == nullptr) {
            result[root->symbol] = code_length;
            return;
        }
        do_tree_traverse(root->left, result, code_length + 1);
        do_tree_traverse(root->right, result, code_length + 1);
    }
};

void generate_fixed_codes() {
    int idx = 0;
    for (unsigned int i = 0b00110000; i <= 0b10111111; i++) {
        fixed_ll_codes.at(idx) = make_pair(8, i);
        idx++;
    }
    for (unsigned int i = 0b110010000; i <= 0b111111111; i++) {
        fixed_ll_codes.at(idx) = make_pair(9, i);
        idx++;
    }
    for (unsigned int i = 0; i <= 0b0010111; i++) {
        fixed_ll_codes.at(idx) = make_pair(7, i);
        idx++;
    }
    for (unsigned int i = 0b11000000; i <= 0b11000111; i++) {
        fixed_ll_codes.at(idx) = make_pair(8, i);
        idx++;
    }
}

vector<pair<u32, u32>> code_lengths_to_code_table(const vector<u32>& lengths) {
    u32 max_len = *max_element(lengths.begin(), lengths.end());
    vector<u32> bl_count(max_len+1, 0);
    //Step 1
    for (u32 i = 0; i < lengths.size(); i++) {
        bl_count[lengths[i]]++;
    }

    //Step 2
    u32 code = 0;
    bl_count[0] = 0;
    vector<u32> next_code(max_len+1, 0);
    for (u32 bits = 1; bits <= max_len; bits++) {
        code = (code + bl_count[bits-1]) << 1;
        next_code[bits] = code;
    }

    //Step 3
    vector<pair<u32, u32>> code_table(lengths.size());
    for (u32 i = 0; i < lengths.size(); i++) {
        u32 len = lengths[i];
        if (len != 0) {
            code_table[i] = make_pair(len, next_code[len]);
            next_code[len]++;
        } else {
            code_table[i] = make_pair(0, 0);
        }
    }
    return code_table;
}

vector<pair<u32, u32>> rle(vector<u32> lengths) {
    const u32 max_run = 138;
    const u32 max_rep = 6;
    size_t n = lengths.size();
    vector<pair<u32, u32>> res;
    for (u32 i = 0; i < n; i++) {
        u32 cur_len = lengths[i];
        if (cur_len != 0) {
            res.push_back(make_pair(lengths[i], 0));
            u32 count = 0;
            while (i < n - 1 && lengths[i+1] == lengths[i] && count < max_rep) {
                count++;
                i++;
            }
            if (count >= 3) {
                res.push_back(make_pair(16, count - 3));
            } else {
                while (count--) {
                    res.push_back(make_pair(cur_len, 0));
                }
            }
        } else {
            u32 count = 1;
            while (i < n - 1 && lengths[i+1] == lengths[i] && count < max_run) {
                count++;
                i++;
            }
            if (count > 10) {
                res.push_back(make_pair(18, count - 11));
            } else if (count >= 3) {
                res.push_back(make_pair(17, count - 3));
            } else {
                while (count--) {
                    res.push_back(make_pair(lengths[i], 0));
                }
            }
        }
    }
    return res;
}

pair<u32, u32> back_ref_search(const array<u8, CHUNK_SIZE>& s, u32 block_size, u32 cur) {
    u32 max_len = 0;
    u32 dist = 0;
    for (int i = cur - 1; i >= 0; i--) {
        if (s[i] == s[cur]) {
            u32 len = 0;
            u32 _dist = cur - i;
            u32 _cur = cur;
            int _i = i;
            while (s[_i] == s[_cur] && len < MAX_MATCH_LENGTH && _cur <= block_size - 1) {
                len++;
                _i++;
                _cur++;
            }
            if (len >= 3 && len > max_len && dist < 32768) {
                max_len = len;
                dist = _dist;
//                return make_pair(len, dist);
            }
//            break;
        }
    }
    return make_pair(max_len, dist);
}

array<u32, 4>* get_LLSymbol(u32 len) {
    for (auto sym: length_code_ranges) {
        if (sym[2] <= len && sym[3] >= len) {
            return new array<u32, 4>(sym);
        }
    }
    return nullptr;
}

array<u32, 4>* get_distSymbol(u32 dist) {
    for (auto sym: dist_code_ranges) {
        if (sym[2] <= dist && sym[3] >= dist) {
            return new array<u32, 4>(sym);
        }
    }
    return nullptr;
}

u32 estimate_block_size(const vector<pair<u32, pair<u32, u32>>>& encoded_stream, const vector<pair<u32, u32>>& ll_codes, const vector<pair<u32, u32>>& dist_codes) {
    u32 size = 0;
    for (u32 i = 0; i < encoded_stream.size(); i++) {
        auto symbol = encoded_stream[i];
        if (symbol.first < 257) {
            // encode single symbol
            auto code = ll_codes[symbol.first];
            size += code.first;
        } else {
            //encode len/dist pair
            auto code = ll_codes[symbol.first];
            size += code.first;
            u32 len = symbol.second.first;
            u32 dist = symbol.second.second;
            array<u32, 4>* lenSym = get_LLSymbol(len);
            array<u32, 4>* distSym = get_distSymbol(dist);
            u32 diff = len - (*lenSym)[2];
            u32 num_offset_bits = (*lenSym)[1];
            if (num_offset_bits > 0) {
                size += num_offset_bits;
            }
            pair<u32, u32> dist_code;
            if (dist_codes.empty()) {
                size += 5;
            } else {
                dist_code = dist_codes[(*distSym)[0]];
            }
            size += dist_code.first;
            diff = dist - (*distSym)[2];
            num_offset_bits = (*distSym)[1];
            if (num_offset_bits > 0) {
                size += num_offset_bits;
            }
        }
    }
    return size;
}

void output_fixed(OutputBitStream& stream, const vector<pair<u32, pair<u32, u32>>>& encoded_stream, bool last_block) {
    stream.push_bit(last_block ? 1 : 0);
    stream.push_bits(1, 2);

    for (u32 i = 0; i < encoded_stream.size(); i++) {
        auto symbol = encoded_stream[i];
        if (symbol.first < 257) {
            auto pair = fixed_ll_codes[symbol.first];
            stream.push_bits(reverse_bits(pair.second, pair.first), pair.first);
        } else {
            //encode len/dist pair
            auto code = fixed_ll_codes[symbol.first];
            stream.push_bits(reverse_bits(code.second, code.first), code.first);
            u32 len = symbol.second.first;
            u32 dist = symbol.second.second;
            array<u32, 4>* lenSym = get_LLSymbol(len);
            array<u32, 4>* distSym = get_distSymbol(dist);
            u32 diff = len - (*lenSym)[2];
            u32 num_offset_bits = (*lenSym)[1];
            if (num_offset_bits > 0) {
                stream.push_bits(diff, num_offset_bits);
            }
            auto dist_code = (*distSym)[0];
            stream.push_bits(reverse_bits(dist_code, 5), 5);
            diff = dist - (*distSym)[2];
            num_offset_bits = (*distSym)[1];
            if (num_offset_bits > 0) {
                stream.push_bits(diff, num_offset_bits);
            }
        }
    }
    stream.push_bits(0, 7);
    if (last_block) {
        stream.flush_to_byte();
    }
}

void output_fixed(OutputBitStream& stream, const vector<u8>& window, const array<u8, CHUNK_SIZE>& contents, u32 block_size, bool last_block) {
    stream.push_bit(last_block ? 1 : 0);
    stream.push_bits(1, 2);


    for (u32 i = 0; i < block_size; i++) {
        if (i >= 2) {
            int pos = i - 2;
            matcher.insert(contents, pos);
        }

        auto len_dist_pair = matcher.find(contents, block_size, i);

        u32 len = len_dist_pair.first;
        u32 dist = len_dist_pair.second;
        auto lenSym = get_LLSymbol(len);
        auto distSym = get_distSymbol(dist);
        if (len > 0 && dist > 0 && lenSym != nullptr && distSym != nullptr) {
            u32 code = (*lenSym)[0];
            auto codeWord = fixed_ll_codes[code];
            stream.push_bits(reverse_bits(codeWord.second, codeWord.first), codeWord.first);
            u32 diff = len - (*lenSym)[2];
            u32 num_offset_bits = (*lenSym)[1];
            if (num_offset_bits > 0) {
                stream.push_bits(diff, num_offset_bits);
            }
            code = (*distSym)[0];
            stream.push_bits(reverse_bits(code, 5), 5);
            diff = dist - (*distSym)[2];
            num_offset_bits = (*distSym)[1];
            if (num_offset_bits > 0) {
                stream.push_bits(diff, num_offset_bits);
            }
            i += (len - 1);
        } else {
            u32 ch = (u32)contents[i];
            auto pair = fixed_ll_codes[ch];
            stream.push_bits(reverse_bits(pair.second, pair.first), pair.first);
        }
    }
    stream.push_bits(0, 7);
    if (last_block) {
        stream.flush_to_byte();
    }
}

void output_dynamic(OutputBitStream& stream, const vector<u8>& history, const array<u8, CHUNK_SIZE>& contents, u32 block_size, bool last_block) {
    vector<u32> ll_freqs(286, 0);
    ll_freqs[256] = 1;
    vector<u32> dist_freqs(30, 0);
    vector<pair<u32, pair<u32, u32>>> encoded_stream;
    for (u32 i = 0;  i < block_size; i++) {
        if (i >= 2) {
            int pos = i - 2;
            matcher.insert(contents, pos);
        }

        if (i > block_size - 3) {
            ll_freqs[contents[i]]++;
            encoded_stream.push_back(make_pair((u32)contents[i], make_pair(0, 0)));
            continue;
        }

        auto len_dist_pair = matcher.find(contents, block_size, i);
        u32 len = len_dist_pair.first;
        u32 dist = len_dist_pair.second;
        array<u32, 4>* lenSym = get_LLSymbol(len);
        array<u32, 4>* distSym = get_distSymbol(dist);
        if (len > 3 && dist > 0 && lenSym != nullptr && distSym != nullptr) {
            ll_freqs[(*lenSym)[0]]++;
            dist_freqs[(*distSym)[0]]++;
            encoded_stream.push_back(make_pair((u32)(*lenSym)[0], len_dist_pair));

            while (--len > 0) {
                matcher.insert(contents, i++);
            }
        } else {
            ll_freqs[contents[i]]++;
            encoded_stream.push_back(make_pair((u32)contents[i], len_dist_pair));
        }
    }

    if (accumulate(dist_freqs.begin(), dist_freqs.end(), 0) <= 0) {
        output_fixed(stream, encoded_stream, last_block);
        return;
    }

    vector<u32> ll_lengths(286, 0);
    for (auto symbol: HuffmanCodingBuilder::package_merge(ll_freqs, 15)) {
        ll_lengths[symbol.first] = symbol.second;
    }
    vector<u32> dist_lengths(30, 0);
    for (auto symbol: HuffmanCodingBuilder::package_merge(dist_freqs, 7)) {
        dist_lengths[symbol.first] = symbol.second;
    }

    vector<pair<u32, u32>> ll_codes = code_lengths_to_code_table(ll_lengths);
    vector<pair<u32, u32>> dist_codes = code_lengths_to_code_table(dist_lengths);

    u32 hlit = 0;
    for (size_t i = ll_lengths.size() - 1; i >= 257; i--) {
        if (ll_lengths[i] != 0) {
            hlit = (u32)(i - 257) + 1;
            break;
        }
    }

    u32 hdist = 0;
    for (size_t i = dist_lengths.size() - 1; i >= 0; i--) {
        if (dist_lengths[i] != 0) {
            hdist = (u32)i;
            break;
        }
    }

    vector<pair<u32, u32>> rle_ll = rle(vector<u32>(ll_lengths.begin(), ll_lengths.begin() + hlit + 257));
    vector<pair<u32, u32>> rle_dist = rle(vector<u32>(dist_lengths.begin(), dist_lengths.begin() + hdist + 1));

    //Freq of cl  symbols
    vector<u32> cl_freqs(19, 0);
    for (auto p: rle_ll) {
        cl_freqs[p.first]++;
    }
    for (auto p: rle_dist) {
        cl_freqs[p.first]++;
    }

    unordered_map<u32, u32> lengths_map = HuffmanCodingBuilder::package_merge(cl_freqs, 7);
    vector<u32> cl_lengths(19,0);
    for (auto p: lengths_map) {
        cl_lengths[p.first] = p.second;
    }
    vector<pair<u32, u32>> cl_codes = code_lengths_to_code_table(cl_lengths);

    vector<u32> reordered_cl_lengths(19, 0);
    for (int i = 0; i < 19; i++) {
        reordered_cl_lengths[i] = cl_lengths[cl_order[i]];
    }

    u32 hclen = 0;
    for (int i = 18; i >= 0; i--) {
        if (reordered_cl_lengths[i] != 0) {
            hclen = (i+1) - 4;
            break;
        }
    }

    // Check size of our compressed block
    vector<pair<u32, u32>> fixed_codes {fixed_ll_codes.begin(), fixed_ll_codes.end()};
    u32 size1 = estimate_block_size(encoded_stream, fixed_codes, {});
    u32 size2 = estimate_block_size(encoded_stream, ll_codes, dist_codes);
    u32 size2_header = 0;
    {
        // until HCLEN
        size2_header += 17;

        size2_header += (hclen + 4) * 3;

        for (int i = 0; i < ll_lengths.size(); i++) {
            size2_header += ll_lengths[i];
        }

        for (int i = 0; i < dist_lengths.size(); i++) {
            size2_header += dist_lengths[i];
        }
    }
    if (size1 < (size2 + size2_header)) {
        output_fixed(stream, encoded_stream, last_block);
        return;
    }

    stream.push_bit(last_block ? 1 : 0);
    stream.push_bits(2, 2);
    stream.push_bits(hlit, 5);
    stream.push_bits(hdist, 5);
    stream.push_bits(hclen, 4);

    //Encode cl code lengths
    for (u32 i = 0; i < hclen + 4; i++) {
        auto code = reordered_cl_lengths[i];
        stream.push_bits(code, 3);
    }

    //Encode ll code lengths
    for (auto p: rle_ll) {
        auto code = cl_codes[p.first];
        //encode symbol prefix code
        stream.push_bits(reverse_bits(code.second, code.first), code.first);
        if (p.first == 18) {
            //encode 7 bits offset
            stream.push_bits(p.second, 7);
        } else if (p.first == 17) {
            //encode 3 bits offset
            stream.push_bits(p.second, 3);
        } else if (p.first == 16) {
            //encode 2 bits offset
            stream.push_bits(p.second, 2);
        }
    }

    //Encode dist code lengths
    for (auto p: rle_dist) {
        auto code = cl_codes[p.first];
        //encode symbol prefix code
        stream.push_bits(reverse_bits(code.second, code.first), code.first);
        if (p.first == 18) {
            //encode 7 bits offset
            stream.push_bits(p.second, 7);
        } else if (p.first == 17) {
            //encode 3 bits offset
            stream.push_bits(p.second, 3);
        } else if (p.first == 16) {
            //encode 2 bits offset
            stream.push_bits(p.second, 2);
        }
    }

    //Encode bitstream
    for (u32 i = 0; i < encoded_stream.size(); i++) {
        auto symbol = encoded_stream[i];
        if (symbol.first < 257) {
            // encode single symbol
            auto code = ll_codes[symbol.first];
            stream.push_bits(reverse_bits(code.second, code.first), code.first);
        } else {
            //encode len/dist pair
            auto code = ll_codes[symbol.first];
            stream.push_bits(reverse_bits(code.second, code.first), code.first);
            u32 len = symbol.second.first;
            u32 dist = symbol.second.second;
            array<u32, 4>* lenSym = get_LLSymbol(len);
            array<u32, 4>* distSym = get_distSymbol(dist);
            u32 diff = len - (*lenSym)[2];
            u32 num_offset_bits = (*lenSym)[1];
            if (num_offset_bits > 0) {
                stream.push_bits(diff, num_offset_bits);
            }
            auto dist_code = dist_codes[(*distSym)[0]];
            stream.push_bits(reverse_bits(dist_code.second, dist_code.first), dist_code.first);
            diff = dist - (*distSym)[2];
            num_offset_bits = (*distSym)[1];
            if (num_offset_bits > 0) {
                stream.push_bits(diff, num_offset_bits);
            }
        }
    }

    auto eob = ll_codes[256];
    stream.push_bits(reverse_bits(eob.second, eob.first), eob.first);
    if (last_block) {
        stream.flush_to_byte();
    }
}

int main(){
    generate_fixed_codes();

    //See output_stream.hpp for a description of the OutputBitStream class
    OutputBitStream stream {cout};

    //Pre-cache the CRC table
    auto crc_table = CRC::CRC_32().MakeTable();

    //Push a basic gzip header
    stream.push_bytes( 0x1f, 0x8b, //Magic Number
        0x08, //Compression (0x08 = DEFLATE)
        0x00, //Flags
        0x00, 0x00, 0x00, 0x00, //MTIME (little endian)
        0x00, //Extra flags
        0x03 //OS (Linux)
    );


    //Note that the types u8, u16 and u32 are defined in the output_stream.hpp header
    array< u8, CHUNK_SIZE > block_contents {};
    u32 block_size {0};
    u32 bytes_read {0};

    char next_byte {};

    u32 crc {};

    if (!cin.get(next_byte)){
        //Empty input?
        stream.push_bit(1);
        stream.push_bits(1, 2);
        stream.push_bits(0, 7);
        stream.flush_to_byte();
    }else{
        bytes_read++;
        //Update the CRC as we read each byte (there are faster ways to do this)
        crc = CRC::Calculate(&next_byte, 1, crc_table); //This call creates the initial CRC value from the first byte read.
        //Read through the input
        while(1){
            block_contents.at(block_size++) = next_byte;

            if (!cin.get(next_byte))
                break;

            bytes_read++;
            crc = CRC::Calculate(&next_byte,1, crc_table, crc); //Add the character we just read to the CRC (even though it is not in a block yet)

            //If we get to this point, we just added a byte to the block AND there is at least one more byte in the input waiting to be written.
            if (block_size == block_contents.size()){
                output_dynamic(stream, history, block_contents, block_size, false);
                block_size = 0;
                matcher.reset();
            }
        }
    }
    //At this point, we've finished reading the input (no new characters remain), and we may have an incomplete block to write.
    if (block_size > 0){
        output_dynamic(stream, history, block_contents, block_size, true);
        block_size = 0;
        matcher.reset();
    }

    //Now close out the bitstream by writing the CRC and the total number of bytes stored.
    stream.push_u32(crc);
    stream.push_u32(bytes_read);

    return 0;
}
