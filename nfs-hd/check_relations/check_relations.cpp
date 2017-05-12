//Usage: ./check_relations <file_relation> <output> <poly.poly> <lpb> <error>

#include <stdio.h>
#include <cassert>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <NTL/ZZXFactoring.h>

using namespace std;
using namespace NTL;

//The polynomial is not irreducible.
#define IRREDUCIBLE "#1"
//The factorization contains not prime factors.
#define PRIME "#2"
//The factorization is false.
#define FACTORIZATION "#3"
//The factorization is not smooth.
#define SMOOTHNESS "#4"

// Read a polynomial in the correct form.
void read_poly(ZZX& a, string poly)
{
  stringstream sline;
  replace(poly.begin(), poly.end(), ',', ' ');
  sline.str(string());
  sline.clear();
  sline << "[" << poly << "]";
  sline >> a;
}

// Convert char of hex to dec.
long hexDec(char c)
{
#ifndef NDEBUG
  assert(c == '-' ||
      c == '0' ||
      c == '1' ||
      c == '2' ||
      c == '3' ||
      c == '4' ||
      c == '5' ||
      c == '6' ||
      c == '7' ||
      c == '8' ||
      c == '9' ||
      c == 'a' ||
      c == 'b' ||
      c == 'c' ||
      c == 'd' ||
      c == 'e' ||
      c == 'f');
#endif // NDEBUG

  if (c == 'a'){
    return 10;
  } else if (c == 'b'){
    return 11;
  } else if (c == 'c') {
    return 12;
  } else if (c == 'd') {
    return 13;
  } else if (c == 'e') {
    return 14;
  } else if (c == 'f') {
    return 15;
  } else if (c == '-') {
    return -1;
  } else {
    return (long) (c - '0');
  }
}

// Convert a hex string in ZZ.
ZZ hexZZ(string hex)
{
  ZZ z;
  z = 0;
  ZZ base;
  base = 1;
  ZZ minus;
  minus = 1;
  ZZ tmp;
  tmp = 0;
  for (size_t i = hex.size(); i-- > 0 ;) {
    tmp = hexDec(hex[i]);
    if (tmp == -1) {
      minus = -1;
    } else {
      z = z + tmp * base;
    }
    base = base * 16;
  }
  return minus * z;
}

void read_poly_hex(ZZX& a, string poly)
{
  char * split;
  string polytmp;
  char * polyc = (char * ) malloc(sizeof(char) * poly.size() + 1);
  copy(poly.begin(), poly.end(), polyc);
  polyc[poly.size()] = '\0';
  split = strtok(polyc, ",");
  long i = 0;
  while (split != NULL) {
    polytmp.assign(split);
    SetCoeff(a, i, hexZZ(polytmp));
    i++;
    split = strtok(NULL, ",");
  }
}

// Read a polynomial inside a relation.
void read_poly_relation(ZZX& a, string line)
{
  string split;
  stringstream sline;
  sline.str(line);
  getline(sline, split, ':');
#ifndef NDEBUG
  assert(split[0] == 'X');
  assert(split[1] == ' ');
  assert(split[2] == '-' || 
      split[2] == '0' ||
      split[2] == '1' ||
      split[2] == '2' ||
      split[2] == '3' ||
      split[2] == '4' ||
      split[2] == '5' ||
      split[2] == '6' ||
      split[2] == '7' ||
      split[2] == '8' ||
      split[2] == '9' ||
      split[2] == 'a' ||
      split[2] == 'b' ||
      split[2] == 'c' ||
      split[2] == 'd' ||
      split[2] == 'e' ||
      split[2] == 'f');
#endif // NDEBUG
  split.erase(0, 2);
  read_poly_hex(a, split);
}

// Read the file *.poly.
void read_poly_file(Vec < ZZX >& f, char * path_poly)
{
  ifstream file_poly (path_poly);
  string line;
  Pair < ZZX, unsigned int > p_tmp;
  Vec< Pair< ZZX, unsigned int > > v_tmp;
  ZZX a;
  while (getline(file_poly, line)) {
    if (line[0] == 'p' && line[1] == 'o' && line[2] == 'l' && line[3] == 'y') {
      p_tmp.b = (int) (line[4] - '0');
      char * split;
      char * linec = (char * ) malloc(sizeof(char) * line.size() + 1);
      copy(line.begin(), line.end(), linec);
      linec[line.size()] = '\0';
      split = strtok(linec, ":");
      split = strtok(NULL, ":");
      line.assign(split);
      line.erase(0, 1);
      read_poly(a, line);
      p_tmp.a = a;
      v_tmp.append(p_tmp);
      free(linec);
    }
  }
  file_poly.close();

  // Stupid sort
  long n = v_tmp.length();
  do {
    long n_tmp = 0;
    for (long i = 1; i < n; i++) {
      if (v_tmp[i - 1].b > v_tmp[i].b) {
        p_tmp = v_tmp[i - 1];
        v_tmp[i - 1] = v_tmp[i];
        v_tmp[i] = p_tmp;
        n_tmp = i;
      }
    }
    n = n_tmp;
  } while (n != 0);

  for (long i = 0; i < v_tmp.length(); i++) {
    f.append(v_tmp[i].a);
  }
}

//Read the lpb.
void read_lpb(Vec <ZZ>& lpb, char * lpbs, long length)
{
  char * split;
  for (long i = 0; i < length; i++) {
    split = strtok(lpbs, ",");
    lpb.append(power2_ZZ(atol(split)));
  }
}

// TODO: Maybe nb_facto not useful.
// Extract factorizations from relation.
void extract_factorisation(Vec < Vec < ZZ > >& factorizations, string line,
    long nb_facto)
{
  string facto;
  Vec < string > factors;
  char * split;
  char * linec = (char * ) malloc(sizeof(char) * line.size() + 1);
  copy(line.begin(), line.end(), linec);
  linec[line.size()] = '\0';
  split = strtok(linec, ":");
  for (long i = 0; i < nb_facto; i++) {
    split = strtok (NULL, ":");
    facto.assign(split);
    factors.append(facto);
  }
  free(linec);

  for (long i = 0; i < factors.length(); i++) {
    Vec < ZZ > v_tmp;
    linec = (char * ) malloc(sizeof(char) * factors[i].size() + 1);
    copy(factors[i].begin(), factors[i].end(), linec);
    linec[factors[i].size()] = '\0';
    split = strtok(linec, ",");
    facto.assign(split);
    v_tmp.append(hexZZ(facto));
    while (split != NULL) {
      split = strtok(NULL, ",");
      if (split != NULL) {
        facto.assign(split);
        v_tmp.append(hexZZ(facto));
      }
    }
    factorizations.append(v_tmp);
    free(linec);
  }
}

/*
 * Return 0 if an element of factor is not prime, 1 if all the elements of
 * factor  are prime.
 */
unsigned int compute_factorization(ZZ& factorization, Vec <ZZ> factor)
{
  factorization = 1;
  unsigned int not_prime = 1;
  //TODO: maybe it is too large.
  long nb_try = 25;
  for (long i = 0; i < factor.length(); i++) {
    factorization *= factor[i];
    if (!ProbPrime(factor[i], nb_try)) {
      not_prime = 0;
    }
  }

  return not_prime;
}

unsigned int is_smooth(Vec <ZZ> factor, ZZ lpb) {
  for (long i = 0; i < factor.length(); i++) {
    if (factor[i] > lpb) {
      return 0;
    }
  }
  return 1;
}

int main(int argc, char ** argv)
{
  ZZX a;

  //Output and error files.
  ofstream filew;
  ofstream fileerr;

  //Input file.
  ifstream filer (argv[1]);

  //Polynomials f read from file_poly.
  Vec <ZZX> f;

  //Large prime bounds.
  Vec <ZZ> lpb;

  if (argc == 6) {
    filew.open(argv[2]);
    read_poly_file(f, argv[3]);
    read_lpb(lpb, argv[4], f.length());
    fileerr.open(argv[5]);
  } else {
    cerr << "Usage error" << endl;
    cerr << "./check_relations <file_relation> <output> <poly.poly> "
      "<lpb> <error>" << endl;
    return 1;
  }

  //Number of relations.
  ZZ nb_rel;
  nb_rel = 0;
  ZZ nb_not_rel;
  nb_not_rel = 0;
  ZZ nb_err_rel;
  nb_err_rel = 0;
  ZZ nb_err_prime;
  nb_err_prime = 0;
  ZZ nb_err_facto;
  nb_err_facto = 0;
  ZZ nb_err_smooth;
  nb_err_smooth = 0;

  string line;

  //Needed to factor a ZZX.
  Vec< Pair< ZZX, long > > factors;
  ZZ c;

  while (getline(filer, line)) {
    if (line[0] != '#') {
      read_poly_relation(a, line);
      assert(content(a) == 1);
      assert(LeadCoeff(a) > 0);
      assert(deg(a) > 0);
      //TODO: if there exists a irreducible test, use it.
      factor(c, factors, a);
      assert(c == 1);

      if (factors.length() == 1) {
          unsigned int error = 0;
          Vec < Vec < ZZ > > factorization;
          extract_factorisation(factorization, line, f.length());
          for (long j = 0; j < f.length(); j++) {
            ZZ norm;
            norm = abs(resultant(a, f[j]));
            ZZ test;
            unsigned int all_prime = compute_factorization(test,
                factorization[j]);
            if (norm != test) {
              fileerr << line << FACTORIZATION << endl;
              error = 1;
              nb_err_rel++;
              nb_err_facto++;
              break;
            }
            if (!all_prime) {
              fileerr << line << PRIME << endl;
              error = 1;
              nb_err_rel++;
              nb_err_prime++;
              break;
            }
            if (!is_smooth(factorization[j], lpb[j])) {
              fileerr << line << SMOOTHNESS << endl;
              error = 1;
              nb_err_rel++;
              nb_err_smooth++;
              break;
            }
          }
          if (!error) {
            filew << line << endl;
            nb_rel++;
          }
      } else {
        fileerr << line << IRREDUCIBLE << endl;
        nb_not_rel++;
      }
    }
  }

  filer.close();
  filew.close();
  if (argc == 6) {
    fileerr.close();
  }

  cout << "Number of relations: " << nb_rel << "." << endl;
  cout << "Number of false relations: " << nb_not_rel << "." <<endl;
  if (nb_err_rel != 0) {
    assert(nb_err_rel == nb_err_prime + nb_err_facto + nb_err_smooth);
    cout << "Number of relations with error: " << nb_err_rel << "." << endl;
    cout << "  contains not prime factor: " << nb_err_prime << "." << endl;
    cout << "  contains bad factorisation: " << nb_err_facto << "." << endl;
    cout << "  is not smooth: " << nb_err_smooth << "." << endl;
  }

  return 0;
}
