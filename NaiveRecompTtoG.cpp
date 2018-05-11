#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <set>
#include <map>
#include <vector>
#include <limits>
#include <chrono>

typedef unsigned int uint;

//
struct SLPVar
{
  SLPVar(uint l, uint r)
  {
    left = l;
    right = r;
  }

  uint left;
  uint right;
};


//
class SLP
{
  std::vector<uint> vars;

public:
  const static int LETTER_OFFSET = 256;

  uint size()
  {
    return vars.size() / 2;
  }

  uint last_id()
  {
    return LETTER_OFFSET + (vars.size()-2) / 2;
  }

  // return id for added variable
  uint push(uint l, uint r)
  {
    
    vars.push_back(l);
    vars.push_back(r);
    return this->last_id();
  }

  bool is_letter(uint id)
  {
    return (id < LETTER_OFFSET);
  }

  uint l_id(uint id);

  uint r_id(uint id);

  int val(uint id, std::string & ret_val);

  void print_val(uint id);

  void print_rules();

  void print_istr(const std::vector<uint> & istr);

  void decomp(std::string & ret_str);

  void decomp(const std::vector<uint> & istr, std::string & ret_str);
};

uint SLP::l_id(uint id)
{
  if (is_letter(id) || id > this->last_id()) {
    return 0;
  }
  return vars[(id-LETTER_OFFSET) * 2];
}

uint SLP::r_id(uint id)
{
  if (is_letter(id) || id > this->last_id()) {
    return 0;
  }
  return vars[(id-LETTER_OFFSET) * 2 + 1];
}

int SLP::val(uint id, std::string & ret_val)
{
  if (is_letter(id)) {
    ret_val = (char)(id);
    return 1;
  } else if (id > this->last_id()) {
    return 0;
  }
  std::string temp_val;
  val(this->l_id(id), temp_val);
  ret_val = temp_val;
  val(this->r_id(id), temp_val);
  ret_val += temp_val;

  return 1;
}

void SLP::print_val(uint id)
{
  std::string val;
  if (this->val(id, val)) {
    std::cout << val;
  } else {
    std::cerr << "error (out of range for id [" << id << "]";
  }
}

void SLP::print_rules()
{
  for (uint id = LETTER_OFFSET; id <= this->last_id(); ++id) {
    std::string val;
    this->val(id, val);
    std::cout << id << " -> " << this->l_id(id) << " " << this->r_id(id) << " : " << val << std::endl;
  }
}

void SLP::print_istr(const std::vector<uint> & istr)
{
#ifdef DEBUG
  std::cout << "size " << istr.size() << std::endl;
#endif
  for (uint i = 0; i < istr.size(); ++i) {
    if (this->is_letter(istr[i])) {
      std::cout << (char)(istr[i]);
    } else {
      std::cout << "[" << istr[i] << "]";
    }
  }
  std::cout << std::endl;
}

void SLP::decomp(std::string & ret_str)
{
  val(this->last_id(), ret_str);
}

void SLP::decomp(const std::vector<uint> & istr, std::string & ret_str)
{
  ret_str = "";
  for (uint i = 0; i < istr.size(); ++i) {
    std::string temp_str;
    val(istr[i], temp_str);
    ret_str += temp_str;
  }
}







uint lsb(uint x)
{
  if (!x) return 0;
  uint k = 1;
  while (!(x & 1)) {
    x = x >> 1;
    ++k;
  }
  return k;
}



void pairComp(std::vector<uint> & istr, SLP & slp)
{
  std::map<uint, int> group;
  std::map< std::pair<uint, uint>, uint > pmap[2];
  group[istr[0]] = 0;
  for (uint i = 1; i < istr.size(); ++i) {
    group[istr[i]] = 0;
    int b;
    std::pair<uint, uint> key;
    if (istr[i-1] < istr[i]) {
      key = std::make_pair(istr[i-1], istr[i]);
      b = 0;
    } else {
      key = std::make_pair(istr[i], istr[i-1]);
      b = 1;
    }
    if (pmap[b].count(key) == 0) {
      pmap[b][key] = 1;
    } else {
      ++(pmap[b][key]);
    }
  }

  // compute the group (0 or 1) on group array
  //    [0, 0, 0]   initially all 0
  // -> [0, -3, -1] process left to right while storing in the later entries the information about which group is better
  // -> [0, 1, 5] put 1 if the value is minus, otherwise 0
  // -> [0, 1, 0]
  std::map< std::pair<uint, uint>, uint >::iterator pmap_it[2];
  pmap_it[0] = pmap[0].begin();
  pmap_it[1] = pmap[1].begin();
  for (std::map<uint, int>::iterator group_it = group.begin(); group_it != group.end(); ++group_it) {
    static const int sign[] = {-1, 1};
    const int g = ((*group_it).second >= 0)? 0 : 1;
    (*group_it).second = g;
    for (int b = 0; b < 2; ++b) {
      while (pmap_it[b] != pmap[b].end() && ((*(pmap_it[b])).first).first == (*group_it).first) {
        group[((*(pmap_it[b])).first).second] += sign[g] * (int)((*(pmap_it[b])).second);
        ++(pmap_it[b]);
      }
    }
  }

  uint count[] = {0, 0};
  for (int b = 0; b < 2; ++b) {
    std::map< std::pair<uint, uint>, uint >::iterator it = pmap[b].begin();
    while ( it != pmap[b].end()) {
      const uint l_id = (b)? ((*it).first).second : ((*it).first).first;
      const uint r_id = (b)? ((*it).first).first : ((*it).first).second;
      const int l_g = group[l_id];
      const int r_g = group[r_id];
      if (l_g == r_g) {
        // delete the element and move to next
        // http://d.hatena.ne.jp/Ox8000ffff/20100225/1267201650
        pmap[b].erase(it++);
      } else {
        // add (*it).second to count[0] if l_g is 0, otherwise add to count[1]
        count[l_g] += (*it).second;
        ++it;
      }
    }
  }
  // determine \Sigma_l to be the set of id's with group 0 if count[0] > count[1], or 1 otherwise.
  const uint sigma_l = (count[0] >= count[1])? 0 : 1;

#ifdef DEBUG
  for (std::map< std::pair<uint, uint>, uint >::iterator it = pmap[0].begin(); it != pmap[0].end(); ++it) {
    std::cout << "pmap[0][<" << ((*it).first).first << ", " << ((*it).first).second << ">] = " << (*it).second << std::endl;
  }
  for (std::map< std::pair<uint, uint>, uint >::iterator it = pmap[1].begin(); it != pmap[1].end(); ++it) {
    std::cout << "pmap[1][<" << ((*it).first).second << ", " << ((*it).first).first << ">] = " << (*it).second << std::endl;
  }
  for (std::map<uint, int>::iterator group_it = group.begin(); group_it != group.end(); ++group_it) {
    std::cout << "group " << (*group_it).first << ", " << (*group_it).second << std::endl;
  }
  std::cout << "count[0], count[1] " << count[0] << ", " << count[1] << std::endl;
  std::cout << "sigma_l " << sigma_l << std::endl;
#endif

  for (int b = 0; b < 2; ++b) {
    std::map< std::pair<uint, uint>, uint >::iterator it = pmap[b].begin();
    while ( it != pmap[b].end()) {
      const uint l_id = (b)? ((*it).first).second : ((*it).first).first;
      const uint r_id = (b)? ((*it).first).first : ((*it).first).second;
      const uint l_g = group[l_id];
      if (l_g != sigma_l) {
        // delete the element and move to next
        // http://d.hatena.ne.jp/Ox8000ffff/20100225/1267201650
        pmap[b].erase(it++);
      } else {
        const uint new_id = slp.push(l_id, r_id);
        (*it).second = new_id;
        ++it;
      }
    }
  }

  uint new_i = 0;
  for (uint i = 1; i <= istr.size(); ++i) {
    if (i == istr.size()) {
      istr[new_i++] = istr[i-1];
      break;
    }
    const uint l_g = group[istr[i-1]];
    const uint r_g = group[istr[i]];
    if (l_g == r_g || l_g != sigma_l) {
      istr[new_i++] = istr[i-1];
      continue;
    } else {
      uint new_id;
      if (istr[i-1] < istr[i]) {
        new_id = pmap[0][std::make_pair(istr[i-1], istr[i])];
      } else {
        new_id = pmap[1][std::make_pair(istr[i], istr[i-1])];
      }
      istr[new_i++] = new_id;
      ++i;
    }
  }
  istr.resize(new_i);
}



void assign_vars_to_blocks(std::map< std::pair<uint, uint>, uint > & bmap, SLP & slp)
{
  for (std::map< std::pair<uint, uint>, uint >::iterator it = bmap.begin(); it != bmap.end(); ) {
    const uint id = ((*it).first).first;
#ifdef DEBUG
    std::cout << "char : " << (char)((*it).first).first << std::endl;
#endif
    uint rep = ((*it).first).second;
    uint maxgap = rep;
    std::map< std::pair<uint, uint>, uint >::iterator scout = it;
    ++scout;
    while (scout != bmap.end() && ((*scout).first).first == id) {
      uint temp_rep = ((*scout).first).second;
      maxgap = std::max(maxgap, temp_rep - rep);
      rep = temp_rep;
      ++scout;
    }

#ifdef DEBUG
    std::cout << "maxgap = " << maxgap << std::endl;
#endif
    const uint beg_id = slp.last_id() + 1;
    uint last_id;
    last_id = slp.push(id, id);
    for (uint remain = maxgap >> 2; remain; remain = remain >> 1) {
      last_id = slp.push(last_id, last_id);
    }

    uint prev_rep = 0;
    uint last_assign = 0;
    for ( ; it != scout; ++it) {
      rep = ((*it).first).second;
      uint remain = rep - prev_rep;
      uint k = lsb(remain);
      if (k == 1) {
        last_id = id;
      } else {
        last_id = beg_id + k - 2;
      }
#ifdef DEBUG
      std::cout << "remain = " << remain << std::endl;
      std::cout << "k = " << k << std::endl;
      std::cout << "last_id = " << last_id << std::endl;
#endif
      for (remain = remain >> k; remain; remain = remain >> 1) {
        ++k;
        if (remain & 1) {
#ifdef DEBUG
          std::cout << "beg_id + k - 2 = " << beg_id + k - 2 << std::endl;
#endif
          last_id = slp.push(last_id, beg_id + k - 2);
        }
      }
      if (last_assign) {
        last_assign = slp.push(last_assign, last_id);
      } else {
        last_assign = last_id;
      }
      bmap[(*it).first] = last_assign;
      prev_rep = rep;
    }
  }
}


void blockComp(std::vector<uint> & istr, SLP & slp)
{
  std::map< std::pair<uint, uint>, uint > bmap;
  uint i = 1;
  uint rep = 1;
  while (true) {
    while (i < istr.size() && istr[i-1] == istr[i]) {
      ++i;
      ++rep;
    }
    if (rep > 1) {
      bmap[std::make_pair(istr[i-1], rep)] = 1;
    }
    if (i >= istr.size()) break;
    ++i;
    rep = 1;
  }

#ifdef DEBUG
  // show block_set
  for (std::map< std::pair<uint, uint>, uint >::iterator it = bmap.begin(); it != bmap.end(); ++it) {
    std::cout << (char)(((*it).first).first) << ":" << ((*it).first).second << std::endl;
  }
#endif

  assign_vars_to_blocks(bmap, slp);
  
  i = 1;
  rep = 1;
  uint new_i = 0;
  while (true) {
    while (i < istr.size() && istr[i-1] == istr[i]) {
      ++i;
      ++rep;
    }
    if (rep > 1) {
#ifdef DEBUG
      std::cout << "bmap = " << bmap[std::make_pair(istr[i-1], rep)] << std::endl;
#endif
      istr[new_i++] = bmap[std::make_pair(istr[i-1], rep)];
    } else {
      istr[new_i++] = istr[i-1];
    }
    if (i >= istr.size()) break;
    ++i;
    rep = 1;
  }
  istr.resize(new_i);
}


int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "arg error" << std::endl;
    return -1;
  }
  std::string filename = argv[1];

  std::ifstream ifs(filename, std::ios::in);
  if (ifs.fail()) {
    std::cerr << "cannot read file " << filename << std::endl;
    return -1;
  }
  std::istreambuf_iterator<char> it(ifs);
  std::istreambuf_iterator<char> last;
  std::string str(it, last);
  ifs >> str;

  std::vector<uint> istr;
  for (uint i = 0; i < str.length(); ++i) {
    istr.push_back((uint)((unsigned char)str[i]));
  }

  const auto startTime = std::chrono::system_clock::now();
  SLP slp;
  while (true) {
    if(istr.size() <= 1) break;
    blockComp(istr, slp);
    if(istr.size() <= 1) break;
    pairComp(istr, slp);
  }
  const auto endTime = std::chrono::system_clock::now();
  const auto timeSpan = endTime - startTime;
  std::cerr << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
  std::cerr << "SLP size: " << slp.size() << std::endl;

  { //// Checking correctness.
    std::cerr << "Checking if SLP can be decompressed correctly..." << std::endl;
    std::string decomp_str;
    slp.decomp(istr, decomp_str);
    // std::cout << decomp_str;
    if (decomp_str == str) {
      std::cerr << "Decompressed correctly." << std::endl;
    } else {
      std::cerr << "Decompression failed." << std::endl;
    }
  }

  return 0;
}

