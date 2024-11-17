// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * cpk_pin.h: 
// * Long-range pseudoknot in human-mRNA: stem1 in 5'UTR, and binding the stem of stem2 from the 3' UTR: 
// *    - stem1: 5' and 3' from 5'UTR;
// *    - stem2: 5' in 5'UTR and 3' from 3'UTR; 
// *
// *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifndef _phylo_h
#define _phylo_h

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <time.h>

using namespace std;

// defines and such  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#define S1_MIN       7     //5, 3
#define S1_MAX       15
#define S2_MIN       6     //5, 3
#define S2_MAX       7
#define L1_MIN       1
#define L1_MAX       2
#define L2_MIN       3
#define L2_MAX       15  //15
#define IS_MIN       0
#define IS_MAX       0
#define MIN_MIS_LEN  6
#define MIN_SUM_LEN  8
#define INIT_BUFF_SIZE		200
#define ENGY_MAX       200
#define CHECK_MAX    30

typedef struct {
  int type;
  int num;
  double engy;
  double delta;
  int start;
  int end;
  int s1;
  int s2;
  int s3;
  int l1;
  int l2;
  int is;
  int mismatch;
  int  mis_s1;
  int  mis_s2;
  int blg;
} pseudoknot_t;

// class definition * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class phyloRna_c {
	// these will hold all the nucleotides that we load in the file
	std::vector<char> list;
	char *file_name;
	char result_name[25];
	char fname[25];		// to store the new file's name
        int cds_start, cds_end;
        int stem_11, stem_12;
	int	stem_21, stem_22;
	int lp_11, lp_12;
	int lp_21, lp_22;
	int is_11, is_12;
        int seq_thrd;
	double max_dist;
	vector <pseudoknot_t> knots;
	int total;
	int allrank;
	FILE *output_file;
	int bulge;
	int mismatch;
	pseudoknot_t knot_engy[ENGY_MAX];
        int verbose_mode; // 0 - simplify; 1 - detail

public:
  static const int S3_MIN = 5;     //stem3-stem
  static const int S3_MAX = 15;     //stem3-stem
  phyloRna_c ();
  phyloRna_c (char *filename, int stem11, int stem12,
	      int stem21, int stem22,
	      int lp11, int lp12, int lp21, int lp22,
	      int is11, int is12);
  ~phyloRna_c ();
  
  char *get_next_line (FILE *fp);
  
  int load_data ();
  int load_result();
  void find_pseudoknot();
  bool base_pair (int p1, int p2);
  void print_knot (pseudoknot_t knot);
  void print_3stems_knot (int pkn, pseudoknot_t knot);
  void timer(int sec, int msec);
  pseudoknot_t extract_result_from_line (char *line);
  pseudoknot_t cal_energy(pseudoknot_t knot);
  bool is_GU_pair(int p1, int p2, int p3, int p4);
  double get_wc_engy(int p1, int p2);
  double get_gu_engy(int p1, int p2, int p3, int p4);
  double get_3d_ends_engy(int p1, int p2, int p3);
  double get_5d_ends_engy(int p1, int p2, int p3);
  int comp_store_engy(pseudoknot_t knot);
  bool comp_to_keep(pseudoknot_t knot);
  void print_engy_knots();
  void print_all_knots();
};


// methods  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


// constructor
phyloRna_c::phyloRna_c () {
	// to be safe, we'll initialize a couple things
	file_name = NULL;
}


// constructor that allows us to set the vlaues for file_name and max_dist
phyloRna_c::phyloRna_c (char *filename, int stem11, int stem12,
			int stem21, int stem22,
			int lp11, int lp12, int lp21, int lp22,
			int is11, int is12) {
	file_name = filename;
	strcpy(fname,file_name);
	strcat(fname,".out");
	stem_11 = stem11;
	stem_12 = stem12;
	stem_21 = stem21;
	stem_22 = stem22;
	lp_11 = lp11;
	lp_12 = lp12;
	lp_21 = lp21;
	lp_22 = lp22;
	is_11 = is11;
	is_12 = is12;
	
	output_file = fopen (fname, "wt");
	  if (output_file == NULL) {
	    printf ("could not open %s for writing\n", fname);
	  }
	bulge = 0;
	mismatch = 0; // percent
}


// destructor
// clear all the lines that were loaded from the pdb file
phyloRna_c::~phyloRna_c () {
	// free all the original line buffers for the amino acids
  fclose (output_file);
}

// load all the sequences from the input file
// assumes all atoms in a sequence will be listed consecutively
int phyloRna_c::load_data (void) {
  string str;
  string acc ="";
  
	// let the user know where we are and what we're doing
	printf ("\n------ load_data (%s) -------\n", file_name);
	if (verbose_mode == 1) {
	  fprintf (output_file, "\n------ load_data (%s) -------\n", file_name);
	}

	FILE *file_pointer;			// pointer the the input file

	// try to open the input file
	file_pointer = fopen (file_name, "rt");

	// bail if we failed to open the file
	if (file_pointer == NULL) {
		printf ("error opening file: %s\n", file_name);
		return -1; }

	char *line;					// the line loaded from the input file
	//get first line - cds
	line = get_next_line (file_pointer);
	//printf("First line: %s \n", line);
	char delimiter = ',';
	for (int i = 0; i < strlen(line); i++) {
	  //cout << str.at(j) << "*";
	  if (line[i] == delimiter) {
	    str = acc;
	    cds_start = atoi(str.c_str());
	    acc = "";
	  } else {
	    acc += line[i];
	  }  
	}
	cds_end = atoi(acc.c_str());

	do {
	  // get the next line from the database file
	  line = get_next_line (file_pointer);
	  
	  // check for EOF
	  if (line == NULL) continue;
	  
	  for (int i = 0; i < strlen(line); i++) {
	    //cout << str.at(j) << "*";
	    list.push_back(line[i]);
	  }
	  
	} while (line != NULL);

	// close the input file
	fclose (file_pointer);

	// tell us about what just happened
	printf ("loaded %ld nucleotides, cds_start %d to cds_end %d\n", list.size(), cds_start, cds_end);
        fprintf (output_file, "loaded %ld nucleotides, cds_start:%d to cds_end:%d\n\n", list.size(), cds_start, cds_end);
	if (verbose_mode == 1) {
	  fprintf (output_file, "loaded %ld nucleotides\n", list.size());
	}

	// no problems encountered

	return list.size();
}



void phyloRna_c::find_pseudoknot (void) {
  int s1, s2, s3, stem3;
	int l1, l2;
	int is;
	int p5_loc, p5_loc_max; // 5 prime
	int p3_loc, p3_loc_max, p3_loc_min; //3 prime
	bool is_1st_knot, is_2nd_knot;
	pseudoknot_t knot;
	bool found, s1_paired, s2_paired;
	int seq_max;
	int A0, A1, A2, A3, A4, A5, A6, A7;
	int A8, A9, A10, A11;
	int i;
	
	total = 0;
	allrank = 0;
	seq_max = list.size();
        is = 0;
	//p5_loc_max to guarantee A5 is in 5'UTR
	p5_loc_max = cds_start - 2*stem_11 - stem_21 - lp_11;
	p5_loc = 0;
	while (p5_loc <= p5_loc_max) {
	  found=false;
	  for (s1 = stem_12;  s1 >= stem_11; s1--) {
            for (s2 = stem_22;  s2 >= stem_21; s2--) {
              for (l1 = lp_11;  l1 <= lp_12; l1++) {

                if ((p5_loc+2*s1+s2+l1) > seq_max) break;
	      

                knot.start = p5_loc;
                knot.s1 = s1;
                knot.s2 = s2;
                knot.l1 = l1;
                knot.l2 = 0;
                knot.is = is;
                //knot.end = p5_loc+2*s1+l1+2*s2+l2+is-1;
                knot.blg = 0;
                knot.mis_s2 = 0;
                
                // find first half of pk
                if (knot.s1+knot.s2<MIN_SUM_LEN) continue;
                A0 = knot.start;
                A1 = A0 + knot.s1-1;
                A2 = A1 + knot.l1+1;
                A3 = A2 + knot.s2-1;
                A4 = A3 + knot.is+1;
                A5 = A4 + knot.s1-1;
                
                if ( !base_pair(A0, A5) ) continue;
                if ( !base_pair(A1, A4) ) continue;
                s1_paired = true;
                for (i=1; i < (knot.s1-1); i++) {
                  if ( !base_pair(A0+i, A5-i) ) {
                    s1_paired = false;
                    break;
                  }
                }

                if (!s1_paired) continue;
                // At this point we get the stem1 of pk
                // The big loop back to pair with stem2 and stem3 is from 3'UTR
                // A7 is between cds_end to list end.
                //fprintf(output_file, "Paired 1st half %d: A0=%d s1=%d s2=%d l1=%d l2=%d A2=%d A3=%d A5=%d\n",
                //        total, A0, s1, s2, l1, l2, A2, A3, A5);
                
                p3_loc_min = cds_end + S3_MIN + knot.s2 +L2_MIN;
                
                for (p3_loc=p3_loc_min; p3_loc <= seq_max; p3_loc++) {
                  found = false;
                  A7 = p3_loc;
                  A6 = A7 - knot.s2 +1;
                  if ( !base_pair(A7, A2) ) continue;
                  if ( !base_pair(A6, A3) ) continue;
                  s2_paired = true;

                  for (i=1; i < (knot.s2-1); i++) {
                    
                    if ( !base_pair(A7-i, A2+i) ) {
                      s2_paired = false;
                      break;               
                    }
                  }
                  if (!s2_paired) continue;
                  
                  // At this point, we got stem2 paired
                  // fprintf(output_file, "Stem2 paired : A2=%d A3=%d s1=%d s2=%d l1=%d l2=%d A6=%d A7=%d\n",
                  //        A2, A3, s1, s2, l1, l2, A6, A7);
                                  
                  A8 = A5 + 1;
                  for (l2 = lp_21;  l2 <= lp_22; l2++) {
                    A11 = A6 - l2 - 1;
                    if (A11 < cds_end) break; // make sure it is on the side of cds_end
                    if ( !base_pair(A8, A11) ) continue;
                    stem3 = 1;
                    for (s3=1; s3<=S3_MAX; s3++) {
                      if ((A8+s3)>cds_start) break;
                      if ((A11-s3)<cds_end) break;
                      if ( !base_pair(A8+s3, A11-s3) ) {
                        //fprintf(output_file, "Stem3 paired %d: A2=%d s1=%d s2=%d s3=%d l1=%d l2=%d A6=%d A7=%d A8=%d A11=%d\n",
                        //         stem3, A2, s1, s2, s3, l1, l2, A6, A7, A8, A11);
                        if (stem3 >= S3_MIN) {found=true;}
                        break;
                      }
                      stem3 +=1;
                    } //for s3
                    if (found) {
                      // we found a pseudoknot with 3 stems
                      total ++;
                      knot.num = total;
                      knot.l2 = l2;
                      knot.s3 = stem3;
                      knot.end = A7;
                      //knots.push_back(knot);
                      //allrank ++;
                      knot = cal_energy(knot);
                      if (knot.engy >18.0) {
                        allrank = comp_store_engy(knot);
                      } 
 
                      found = false;
                    } // if found
                    //continue with different l2 to find any more of pk
                  } // for l2
		  //continue with different p3 (A7) to find any more of pk
                } // for p3_loc
                //new round of l1 to search for first half of pk
              } //l1
	      //new round of l1 to search for first half of pk
            } // s2
	  } // s1
	  p5_loc ++;
	} // while
}


bool phyloRna_c::base_pair (int p1, int p2) {

	switch (list[p1]) {
	  case 'a':
	  case 'A':
			if (list[p2]=='t' || list[p2]=='T') 
				return true;
			else if (list[p2]=='u' || list[p2]=='U') 
				return true;
			else 
				return false;
			break;
	  case 'u':
	  case 'U':
	  case 't':
	  case 'T':
			if (list[p2]=='a' || list[p2]=='A') 
			  return true;
			else if (list[p2]=='g' || list[p2]=='G') 
			  return true;
			else 
			  return false;
			break;		
	  case 'g':
	  case 'G':
			if (list[p2]=='c' || list[p2]=='C') 
			  return true;
			else if (list[p2]=='u' || list[p2]=='U') 
			  return true;
			else if (list[p2]=='t' || list[p2]=='T') 
			  return true;
			else 
			  return false;
			break;		
	  case 'c':
	  case 'C':
			if (list[p2]=='g' || list[p2]=='G') return true;
			else return false;
			break;	
	  default:
		  return false;
		  break;
	} // switch
}


pseudoknot_t phyloRna_c::cal_energy(pseudoknot_t knot) {
	int i;
	double energy;
	int A0, A1, A2, A3, A4, A5, A6, A7;
        int A8, A9, A10, A11;
	
	A0 = knot.start;
	A1 = A0 + knot.s1-1;
	A2 = A1 + knot.l1+1;
	A3 = A2 + knot.s2-1;
	A4 = A3 + knot.is+1;
	A5 = A4 + knot.s1-1;

        A7 = knot.end;
	A6 = A7 - knot.s2 +1;

        A8 = A5 + 1;
        A9 = A8 + knot.s3 -1;
        
        A11 = A6 - knot.l2 - 1;
        A10 = A11 - knot.s3 +1;
        
	energy = 0.0;
	// Calculate 5' dangling ends
	if (A0>1)
		energy += get_5d_ends_engy(A5, A0, A0-1);
	//fprintf(output_file, "rank=%d 5d engy=%f \n", knot.num, get_5d_ends_engy(A5, A0, A0-1));
	// Calculate 3' dangling ends
	energy += get_3d_ends_engy(A7, A2, A7+1);
	//fprintf(output_file, "rank=%d 3d  engy=%f \n", knot.num,get_3d_ends_engy(A7, A2, A7+1));
	for (i=0; i < (knot.s1-1); i++) {
	  if ( is_GU_pair(A0+i, A5-i, A0+i+1, A5-i-1) ) {
	    energy += get_gu_engy(A0+i, A5-i, A0+i+1, A5-i-1);
	    //fprintf(output_file, "rank=%d engy=%f s1 GU %c-%c:%c-%c\n", 
	    //	knot.num,get_gu_engy(A0+i, A5-i, A0+i+1, A5-i-1), list[A0+i],list[A5-i],list[A0+i+1],list[A5-i-1]);
		} else {
	    energy += get_wc_engy(A0+i,A0+i+1);
	    //fprintf(output_file, "rank=%d engy=%f s1 %c:%c\n", 
	    //	knot.num,get_wc_engy(A0+i,A0+i+1), list[A0+i],list[A0+i+1]);
		} // if
	} //  for

	for (i=0; i < (knot.s2-1); i++) {
	  if ( is_GU_pair(A6+i, A3-i, A6+i+1, A3-i-1) ) {
	    energy += get_gu_engy(A6+i, A3-i, A6+i+1, A3-i-1);
	    //fprintf(output_file, "rank=%d engy=%f s2 GU %c-%c:%c-%c\n", 
	    //	knot.num,get_gu_engy(A6+i, A3-i, A6+i+1, A3-i-1), list[A6+i],list[A3-i],list[A6+i+1],list[A3-i-1]);
	  } else {
	    energy += get_wc_engy(A6+i, A6+i+1);
	    //fprintf(output_file, "rank=%d engy=%f s2 %c:%c\n", 
	    //	knot.num,get_wc_engy(A6+i, A6+i+1),list[A6+i],list[A6+i+1]);
	  } // if
	} // for

        
	if (knot.is == 0) {
	  if ( is_GU_pair(A1, A4, A6, A3) ) {
	    energy += 0.5*get_gu_engy(A1, A4, A6, A3);
	    //fprintf(output_file, "rank=%d engy=%f 1/2 %c-%c:%c-%c\n", 
	    //	knot.num,0.5*get_gu_engy(A1, A4, A6, A3),list[A1],list[A4],list[A6],list[A3]);
	  } else {
	    energy += 0.5*get_wc_engy(A1, A6);
	    //fprintf(output_file, "rank=%d engy=%f 1/2 %c:%c\n", 
	    //	knot.num,0.5*get_wc_engy(A1, A6),list[A1],list[A6]);
	  } // if	  
	}
	
	knot.engy = energy;
	return knot;
}

int phyloRna_c::comp_store_engy(pseudoknot_t knot) {
  int i, j;
  int lp;
  int option;
  int pos;
  
  // if (allrank < ENGY_MAX) 
  //   lp = allrank;
  // else
  //   lp = ENGY_MAX;

  option = 0;
  if (allrank !=0){
    // allrank keep how much we have in the kont_engy array
    if (allrank < ENGY_MAX) 
      lp = allrank;
    else
      lp = ENGY_MAX;
    // search the list to find out any pk have the same start and same end
    // compare to keep the lowest engy  
    for (i=0; i<lp; i++) {
      if ((knot.start == knot_engy[i].start) && (abs(knot.end - knot_engy[i].end)<=10)) {
        if (knot.engy <= knot_engy[i].engy) {
          // do nothing and return
          option = 1;
        } else {
          // remove the element and shift the rest of the pks
          option = 2;
          pos = i;
        }
      }
    }
  }
  
  if (option != 1) {
    if (option == 2) {
      // shift value start from knot_engy[pos] - copy the next one over
      for (i=pos; i<lp-1; i++){ // now the last one [lp-2] has the value of [lp-1]
        knot_engy[i] = knot_engy[i+1];
      }
      // allrank eventually unchange
      allrank --;
    }

    allrank ++;
    // allrank keep how much we have in the kont_engy array
    if (allrank < ENGY_MAX) 
      lp = allrank;
    else
      lp = ENGY_MAX;
    
    knot_engy[lp-1].engy = 0;
    for (i=0; i<lp; i++) {
      if (knot.engy > knot_engy[i].engy) {
        for (j=lp-1; j>i; j--) {
          if ( j-1>=0 ) {
            knot_engy[j] = knot_engy[j-1];
          }
      } //for 
        knot_engy[i] = knot;
        break;
      } // if
    } // for
  }

  return allrank;      
}


bool phyloRna_c::comp_to_keep(pseudoknot_t knot) {
  int i, j;
  pseudoknot_t tmp;
  int lp1, lp2;

  lp2=knots.size();
 
  if (lp2<CHECK_MAX)
    lp1 = 0;
  else
    lp1 = lp2 - CHECK_MAX;


  for (i=lp1; i<lp2; i++) {
    if (knot.end > knots[i].end)
      continue;
    else {
      if (knot.start >= knots[i].start) {
	if (knot.is != knots[i].is)
	  continue;
	else {
	  if (((knot.s1+knot.l2)==(knots[i].s1+knots[i].l2)) &&
	      ((knot.s2+knot.l1)==(knots[i].s2+knots[i].l1)))
	    return false;
	}
      } 
    } // else
  } // for
  return true;
}


void phyloRna_c::print_engy_knots() {
  int i, j;
  int lp;
  int rank;
  bool keep;
  int A2, A5, A6, B2, B5, B6;
	

  //printf("knot_num = %d \n", knot_num);
  //fprintf(output_file,"\n\n-------------------- Energy Rank ------------------------------\n\n");
  if (allrank < ENGY_MAX) // rank only 100 pk 
    lp = allrank;
  else
    lp = ENGY_MAX;

  rank = 1;

  //  fprintf(output_file, "\n Collect %d pseudoknots! \n", lp);

  for (i=0; i<lp; i++) {
    /* fprintf(output_file, "\nEnergy %d: -%f: Start=%d S1=%d S2=%d L1=%d L2=%d L3=%d End=%d\n\n", */

    /* i, knot_engy[i].engy, knot_engy[i].start+1, knot_engy[i].s1, knot_engy[i].s2, knot_engy[i].l1, knot_engy[i].l2, knot_engy[i].is, knot_engy[i].end+1); */
    /* fprintf(output_file, "\n # %d Energy= -%f \n", i, knot_engy[i].engy); */
    keep = true;
    if (knot_engy[i].engy < 18.0) continue;
    for (j=0; j<i; j++) {
       if ((knot_engy[i].start >= knot_engy[j].start)&&(knot_engy[i].end <= knot_engy[j].end)&&(knot_engy[i].engy <= knot_engy[j].engy)) {
	 A2 = knot_engy[i].start + knot_engy[i].s1 +knot_engy[i].l1;
	 A5 = A2 + knot_engy[i].s2+knot_engy[i].is+knot_engy[i].s1-1;
	 A6 = A5 + knot_engy[i].l2+1;
	 B2 = knot_engy[j].start + knot_engy[j].s1 +knot_engy[j].l1;
	 B5 = B2 + knot_engy[j].s2+knot_engy[j].is+knot_engy[j].s1-1;
	 B6 = B5 + knot_engy[j].l2+1;
	 if ((B2<=A2)&&(A5<=B5)&&(A6>=B6)) {
           //fprintf(output_file, "\nEnergy = %f  start= %d A2,A5,B2,B5: %d %d, %d %d s1=%d s1=%d \n", 
	   //	   knot_engy[i].engy, knot_engy[i].start, A2, A5, B2, B5,knot_engy[i].s1,knot_engy[j].s1 );
	   keep = false;
	   break;
	 }
       }
    }
    if (!keep) continue;
    //	fprintf(output_file, "Energy %d: -%f: org_rank=%d start=%d s1=%d s2=%d l1=%d l2=%d is=%d end=%d\n",
    //		i+1, knot_engy[i].engy, knot_engy[i].num, knot_engy[i].start+1, knot_engy[i].s1, knot_engy[i].s2, knot_engy[i].l1, knot_engy[i].l2, knot_engy[i].is, knot_engy[i].end+1);
    fprintf(output_file, "\nEnergy %d: -%f: Start=%d S1=%d S2=%d L1=%d L2=%d L3=%d End=%d\n\n",

    rank, knot_engy[i].engy, knot_engy[i].start+1, knot_engy[i].s1, knot_engy[i].s2, knot_engy[i].l1, knot_engy[i].l2, knot_engy[i].is, knot_engy[i].end+1);
    print_knot(knot_engy[i]);
    rank += 1;
  }
}



bool phyloRna_c::is_GU_pair(int p1, int p2, int p3, int p4) { // p1-p2:p3-p4 (p3=p1+1, p4=p2-1)

  bool value=false;

	switch (list[p1]) {
	  case 'a':
	  case 'A':
	    switch (list[p3]) { //p1+1
		      case 'a':
		      case 'A':
			        value = false;
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			   if ((list[p4]=='g')||(list[p4]=='G'))
			        value = true;
			   else 
			        value = false;
			   
			   break;		
		      case 'g':
		      case 'G':
			   if ((list[p4]=='u')||(list[p4]=='U')||(list[p4]=='t')||(list[p4]=='T'))
			        value = true;
			   else 
			        value = false;
			   
		                break;		
		      case 'c':
		      case 'C':
			        value = false;		
		                  break;	
		      default:
		                  break;
		   } // switch

		   break;
	  case 'u':
	  case 'U':
	  case 't':
	  case 'T':
		  if ((list[p2] == 'a') || (list[p2] == 'A')) {
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value = false;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			   if ((list[p4]=='g')||(list[p4]=='G'))
			        value = true;
			   else 
			        value = false;
			   
		           break;		
		      case 'g':
		      case 'G':
			   if ((list[p4]=='u')||(list[p4]=='U')||(list[p4]=='t')||(list[p4]=='T'))
			        value = true;
			   else 
			        value = false;
			   
		           break;		
		      case 'c':
		      case 'C':
			        value = false;		
		                  break;	
		      default:
		                  break;
		   } // switch
		  } else if ((list[p2] == 'g') || (list[p2] == 'G')){
		    value = true;
		  }
		  break;		
	  case 'g':
	  case 'G':
		  if ((list[p2] == 'c') || (list[p2] == 'C')){
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value = false;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			   if ((list[p4]=='g')||(list[p4]=='G'))
			        value = true;
			   else 
			        value = false;
		                break;		
		      case 'g':
		      case 'G':
			   if ((list[p4]=='u')||(list[p4]=='U')||(list[p4]=='t')||(list[p4]=='T'))
			        value = true;
			   else 
			        value = false;
		                break;		
		      case 'c':
		      case 'C':
			        value = false;		
		                  break;	
		      default:
		                  break;
		   } // switch
		  } else if ((list[p2] == 'u') || (list[p2] == 'U')||(list[p2]=='t')||(list[p2]=='T')){
		    value = true;
		  }
			break;		
	  case 'c':
	  case 'C':
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value = false;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			   if ((list[p4]=='g')||(list[p4]=='G'))
			        value = true;
			   else 
			        value = false;
		                break;		
		      case 'g':
		      case 'G':
			   if ((list[p4]=='u')||(list[p4]=='U')||(list[p4]=='t')||(list[p4]=='T'))
			        value = true;
			   else 
			        value = false;
		                break;		
		      case 'c':
		      case 'C':
			        value = false;		
		                  break;	
		      default:
		                  break;
		   } // switch
		
			break;	
	  default:
		  break;
	} // switch
	return value;
}

double phyloRna_c::get_wc_engy(int p1, int p2) { 
  double value;

        value = 0;
	switch (list[p1]) {
	  case 'a':
	  case 'A':
	           switch (list[p2]) {
		      case 'a':
		      case 'A':
			        value=0.93;
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=1.10;		
		                break;		
		      case 'g':
		      case 'G':
			        value=2.08;		
		                break;		
		      case 'c':
		      case 'C':
			        value=2.24;		
		                  break;	
		      default:
		                  break;
		   } // switch

		   break;
	  case 'u':
	  case 'U':
	  case 't':
	  case 'T':

	           switch (list[p2]) {
		      case 'a':
		      case 'A':
			        value=1.33;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.93;		
		                break;		
		      case 'g':
		      case 'G':
			        value=2.11;		
		                break;		
		      case 'c':
		      case 'C':
			        value=2.35;		
		                  break;	
		      default:
		                  break;
		   } // switch

		
			break;		
	  case 'g':
	  case 'G':
	           switch (list[p2]) {
		      case 'a':
		      case 'A':
			        value=2.35;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=2.24;		
		                break;		
		      case 'g':
		      case 'G':
			        value=3.26;		
		                break;		
		      case 'c':
		      case 'C':
			        value=3.42;		
		                  break;	
		      default:
		                  break;
		   } // switch
		
			break;		
	  case 'c':
	  case 'C':
	           switch (list[p2]) {
		      case 'a':
		      case 'A':
			        value=2.11;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=2.08;		
		                break;		
		      case 'g':
		      case 'G':
			        value=2.36;		
		                break;		
		      case 'c':
		      case 'C':
			        value=3.26;		
		                  break;	
		      default:
		                  break;
		   } // switch
		
			break;	
	  default:
		  ;
		  break;
	} // switch
	return value;
}

double phyloRna_c::get_gu_engy(int p1, int p2, int p3, int p4) { // p1-p2:p3-p4 p2=p1+1, p4=p2-1
  double value;

        value = 0;
	switch (list[p1]) {
	  case 'a':
	  case 'A':
	           switch (list[p3]) {
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=1.36;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.55;		
		                break;			
		      default:
		                  break;
		   } // switch

		   break;
	  case 'u':
	  case 'U':
	  case 't':
	  case 'T':
	    if ((list[p2] == 'a') || (list[p2] == 'A')) {
	           switch (list[p3]) {
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=1.27;		
		                break;		
		      case 'g':
		      case 'G':
			        value=1.0;		
		                break;			
		      default:
		                  break;
		   } // switch
	    } else if ((list[p2] == 'g') || (list[p2] == 'G')) {
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=1.0;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			if ((list[p4] == 'a')||(list[p4] == 'A')) {
			        value=0.55;
			} else if ((list[p4] == 'g')||(list[p4] == 'G')) {
			        value=-0.47;
			}
		                break;		
		      case 'g':
		      case 'G':
			if ((list[p4] == 'u')||(list[p4] == 'U')||(list[p4] == 't')||(list[p4] == 'T')) {
			        value=-0.3;
			} else if ((list[p4] == 'c')||(list[p4] == 'C')) {
			        value=1.41;
			}			        		
			break;		
		      case 'c':
		      case 'C':
			        value=1.53;		
		                  break;	
		      default:
		                  break;
		   } // switch
	    }
	    break;		
	  case 'g':
	  case 'G':
	    if ((list[p2] == 'c') || (list[p2] == 'C')) {
	           switch (list[p3]) {
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=2.51;		
		                break;		
		      case 'g':
		      case 'G':
			        value=1.53;		
		                break;			
		      default:
		                  break;
		   } // switch
	    } else if ((list[p2] == 'u') || (list[p2] == 'U')||(list[p2] == 't')||(list[p2] == 'T')) {
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=1.27;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			if ((list[p4] == 'a')||(list[p4] == 'A')) {
			        value=1.36;
			} else if ((list[p4] == 'g')||(list[p4] == 'G')) {
			        value=-1.29;
			}
		                break;		
		      case 'g':
		      case 'G':
			if ((list[p4] == 'u')||(list[p4] == 'U')||(list[p4] == 't')||(list[p4] == 'T')) {
			        value=-0.47;
			} else if ((list[p4] == 'c')||(list[p4] == 'C')) {
			        value=2.11;
			}			        		
			break;		
		      case 'c':
		      case 'C':
			        value=2.51;		
		                  break;	
		      default:
		                  break;
		   } // switch
	    }
	    break;		
	  case 'c':
	  case 'C':
	           switch (list[p3]) {
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=2.11;		
		                break;		
		      case 'g':
		      case 'G':
			        value=1.41;		
		                break;			
		      default:
		                  break;
		   } // switch
		
			break;	
	  default:
		  ;
		  break;
	} // switch
	return value;
}


double phyloRna_c::get_3d_ends_engy(int p1, int p2, int p3) {
  double value;

  value = 0;
	switch (list[p1]) {
	  case 'a':
	  case 'A':
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=0.8;
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.6;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.8;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.5;		
		                  break;	
		      default:
		                  break;
		   } // switch

		   break;
	  case 'u':
	  case 'U':
	  case 't':
	  case 'T':
		  if ((list[p2] == 'a') || (list[p2] == 'A')) {
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=0.7;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.1;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.7;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.1;		
		                  break;	
		      default:
		                  break;
		   } // switch
		  } else if ((list[p2] == 'g') || (list[p2] == 'G')){
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=0.7;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.1;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.7;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.1;		
		                  break;	
		      default:
		                  break;
		   }
		  }
		
			break;		
	  case 'g':
	  case 'G':
		  if ((list[p2] == 'c') || (list[p2] == 'C')){
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=1.1;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.6;		
		                break;		
		      case 'g':
		      case 'G':
			        value=1.3;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.4;		
		                  break;	
		      default:
		                  break;
		   } // switch
		  } else if ((list[p2] == 'u') || (list[p2] == 'U')){
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=0.8;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.6;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.8;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.5;		
		                  break;	
		      default:
		                  break;
		   } // switch
		  }
			break;		
	  case 'c':
	  case 'C':
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=1.7;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=1.2;		
		                break;		
		      case 'g':
		      case 'G':
			        value=1.7;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.8;		
		                  break;	
		      default:
		                  break;
		   } // switch
		
			break;	
	  default:
		  break;
	} // switch
	return value;
     
}

double phyloRna_c::get_5d_ends_engy(int p1, int p2, int p3) {
  double value;

  value = 0;
	switch (list[p1]) {
	  case 'a':
	  case 'A':
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=0.3;
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.2;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.2;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.1;		
		                  break;	
		      default:
		                  break;
		   } // switch

		   break;
	  case 'u':
	  case 'U':
	  case 't':
	  case 'T':
		  if ((list[p2] == 'a') || (list[p2] == 'A')) {
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=0.3;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.2;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.4;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.3;		
		                  break;	
		      default:
		                  break;
		   } // switch
		  } else if ((list[p2] == 'g') || (list[p2] == 'G')){
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=0.3;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.2;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.4;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.3;		
		                  break;	
		      default:
		                  break;
		   }
		  }		
			break;		
	  case 'g':
	  case 'G':
		  if ((list[p2] == 'c') || (list[p2] == 'C')){
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=0.5;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.1;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.2;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.3;		
		                  break;	
		      default:
		                  break;
		   } // switch
		  } else if ((list[p2] == 'u') || (list[p2] == 'U')){
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=0.3;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.2;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.2;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.1;		
		                  break;	
		      default:
		                  break;
		   } // switch
		  }
			break;		
	  case 'c':
	  case 'C':
	           switch (list[p3]) {
		      case 'a':
		      case 'A':
			        value=0.2;		
			        break;
		      case 'u':
		      case 'U':
		      case 't':
		      case 'T':
			        value=0.0;		
		                break;		
		      case 'g':
		      case 'G':
			        value=0.0;		
		                break;		
		      case 'c':
		      case 'C':
			        value=0.3;		
		                  break;	
		      default:
		                  break;
		   } // switch
		   break;	
	  default:
		  break;
	} // switch
	return value;
}


// bool phyloRna_c::is_first_half_pk (pseudoknot_t &knot) {
// 	int A0, A1, A2, A3, A4, A5, A6, A7;
// 	int i;

// 	// when we have s1 s2 and l1
// 	if (knot.s1+knot.s2<MIN_SUM_LEN) return false;
// 	A0 = knot.start;
// 	A1 = A0 + knot.s1-1;
// 	A2 = A1 + knot.l1+1;
// 	A3 = A2 + knot.s2-1;
// 	A4 = A3 + knot.is+1;
// 	A5 = A4 + knot.s1-1;
	
// 	int mis_s1=0;
// 	int mis_s2=0;
	

// 	if ( !base_pair(A0, A5) ) return false;
// 	if ( !base_pair(A1, A4) ) return false;

// 	for (i=1; i < (knot.s1-1); i++) {
// 		if ( !base_pair(A0+i, A5-i) ) {
// 		  if ((mismatch!=0)&&(knot.s1>MIN_MIS_LEN)) { // only allow mismatch when stem >5
// 				mis_s1 ++;
// 				if (mis_s1 > mismatch) return false;
// 			} else { // no bulge or mismatch allowed
// 				return false;
// 			}
// 		}
// 	}

	// for (i=1; i < (knot.s2-1); i++) {
	// 	if ( !base_pair(A6+i, A3-i) ) {
	// 	  if ((mismatch!=0)&&(knot.s2>MIN_MIS_LEN)) { // only allow mismatch when stem >5
	// 			mis_s2 ++;
	// 			if ((mis_s1+mis_s2) > mismatch) return false;
	// 		} else { // no bulge or mismatch allowed
	// 			return false;
	// 		}
	// 	}
	// }
			
// 	knot.mismatch =  mis_s1+mis_s2;
// 	knot.mis_s1 = mis_s1;
// 	knot.mis_s2 = mis_s2;
	
// 	return true;
// }



void phyloRna_c::print_all_knots() {
  int i, j;
  int lp;
  int rank;
  //int arr_size;

  if (allrank < ENGY_MAX) // rank only 100 pk 
    lp = allrank;
  else
    lp = ENGY_MAX;

  rank = 1;
  
  //  arr_size = sizeof(knot_engy) / sizeof(knot_engy[0]);
  printf("Printing %d pseudoknots.\n", lp);
  for (i=0; i<lp; i++) {
    print_3stems_knot(i, knot_engy[i]);
  }
}

void phyloRna_c::print_3stems_knot (int pkn, pseudoknot_t knot) {
   int A0, A1, A2, A3, A4, A5, A6, A7;
   int A8, A9, A10, A11, l3;
   int i, seq_max;
   int pcase;

        seq_max = list.size();
	A0 = knot.start;
	A1 = A0 + knot.s1-1;
	A2 = A1 + knot.l1+1;
	A3 = A2 + knot.s2-1;
	A4 = A3 + knot.is+1;
	A5 = A4 + knot.s1-1;

        A7 = knot.end;
	A6 = A7 - knot.s2+1;
        A8 = A5 + 1;
	A9 = A8 + knot.s3-1;
        A11 = A6 -knot.l2 - 1;
        A10 = A11 - knot.s3 +1;

        pcase = 0;

        fprintf(output_file, "PK# %d: Engy=%f start=%d s1=%d s2=%d s3=%d l1=%d l2=%d end=%d\n",
            pkn+1, knot.engy, knot.start, knot.s1, knot.s2, knot.s3,
            knot.l1, knot.l2, knot.end);
        fprintf(output_file, "         A0=%d A1=%d A2=%d A3=%d A5=%d A6=%d A7=%d(%d) A9=%d A10=%d(%d) A11=%d \n\n",
                A0, A1, A2, A3, A5, A6, A7, A7-cds_end+1, A9, A10, A10-cds_end+1, A11);        
	// First line: L2, L1
	//printf ("                    ");
	//	fprintf (output_file, "Bulge is %d, Pcase is %d\n", knot.blg, pcase);
	fprintf (output_file, "         ");
        fprintf (output_file, "                    ");

	if (knot.l2 >=100) {
	  fprintf (output_file, "%d nt", knot.l2);
	} else {
	  for (i=1; i<=knot.l2; i++) {
	    //printf("%c", list[A5+i]);
	      fprintf (output_file, "%c", toupper(list[A11+i]));
	  } //for
	} // else

       //printf("       ")
        fprintf (output_file, "       ");
        
	for (i=1; i<=knot.l1; i++) {
	 //printf("%c", list[A1+i]);
          fprintf (output_file, "%c", toupper(list[A1+i]));
	}

       //printf("\n\n");
        fprintf (output_file, "\n\n");

	// Second line: Stem3, Stem1-R, IS, Stem2-R
	//printf ("                    ");
	fprintf (output_file, "                    ");

	for (i=0; i<knot.s3; i++) {
          fprintf (output_file, "%c", toupper(list[A9-i]));
	}
        fprintf (output_file, "       ");

	for (i=0; i<knot.s1; i++) {
	 //printf("%c", list[A5-i]);
          fprintf (output_file, "%c", toupper(list[A5-i]));
	}

       //printf("  ");
        fprintf (output_file, "  ");

	for (i=1; i<=knot.is; i++) {
	 //printf("%c", list[A4-i]);
          fprintf (output_file, "%c", toupper(list[A4-i]));
	}

       //printf("  ");
        fprintf (output_file, "  ");

	for (i=0; i<knot.s2; i++) {
	 //printf("%c", list[A3-i]);
          fprintf (output_file, "%c", toupper(list[A3-i]));
	}

       //printf("\n\n");
        //fprintf (output_file, "\n\n");
        l3 = A10 - A9 -1;
        fprintf (output_file, "\n");
        fprintf (output_file, "       %d nts", l3);
        fprintf (output_file, "\n");

	// Third line: Stem3, *Stem1, Stem1-L, Stem2-L
        fprintf (output_file, "                    ");
                
	for (i=0; i<knot.s3; i++) {
	 //printf("%c", list[A0+i]);
          fprintf (output_file, "%c", toupper(list[A10+i]));
	}

        fprintf (output_file, "     ");

        fprintf (output_file, "**");

	for (i=0; i<knot.s1; i++) {
	 //printf("%c", list[A0+i]);
          fprintf (output_file, "%c", toupper(list[A0+i]));
	}

       //printf("  ");
        fprintf (output_file, "  ");

	for (i=1; i<=knot.is; i++) {
	 //printf(" ");
	  fprintf (output_file, " ");
	}

       //printf("  ");
        fprintf (output_file, "  ");
	
	for (i=0; i<knot.s2; i++) {
          fprintf (output_file, "%c", toupper(list[A6+i]));
	}

	// 2 more after A7
        fprintf (output_file, "*");
	for (i=1; i<3; i++) {
          if ((A7+i)<seq_max) 
            fprintf (output_file, "%c", toupper(list[A7+i]));
	}

       //printf("\n\n");	
        fprintf (output_file, "\n\n");	
}

void phyloRna_c::print_knot (pseudoknot_t knot) {
   int A0, A1, A2, A3, A4, A5, A6, A7;
   int i;
   int pcase;

	A0 = knot.start;
	A1 = A0 + knot.s1-1;
	A2 = A1 + knot.l1+1;
	A3 = A2 + knot.s2-1;
	A4 = A3 + knot.is+1;
	A5 = A4 + knot.s1-1;
	A6 = A5 + knot.l2+1;
	A7 = A6 + knot.s2-1;

	if (knot.blg == 0) 
	  pcase = 0;
	else {
	  if ((knot.blg >= A0) &&(knot.blg <= A1)) {
	    pcase = 1;
	  } else  if ((knot.blg >= A2) &&(knot.blg <= A3)) {
	    pcase = 2;
	  } else if ((knot.blg >= A4) &&(knot.blg <= A5)) {
	    pcase = 3;
	  } else if ((knot.blg >= A6) &&(knot.blg <= A7)) {
	    pcase = 4;	
	  }
	}

	// First line: L2, L1
	//printf ("                    ");
	//	fprintf (output_file, "Bulge is %d, Pcase is %d\n", knot.blg, pcase);
	fprintf (output_file, "                    ");
	if (knot.l2 >=100) {
	  fprintf (output_file, "%d nt", knot.l2);
	} else {
	  for (i=1; i<=knot.l2; i++) {
	    //printf("%c", list[A5+i]);
	    switch (pcase) {
	    case 0:
	    case 1:
	    case 2:
	    case 3:
	      fprintf (output_file, "%c", toupper(list[A5+i+bulge]));
	      break;
	    case 4:
	      fprintf (output_file, "%c", toupper(list[A5+i]));
	      break;
	    } // switch
	  } //for
	} // else

       //printf("       ");
        fprintf (output_file, "       ");
	
	for (i=1; i<=knot.l1; i++) {
	 //printf("%c", list[A1+i]);
	  switch (pcase) {
	    case 0:
	    case 2:
	    case 3:
	    case 4:
	      fprintf (output_file, "%c", toupper(list[A1+i]));
	      break;
	    case 1:
	      fprintf (output_file, "%c", toupper(list[A1+i+bulge]));
	      break;
	  } // switch
	}

       //printf("\n\n");
        fprintf (output_file, "\n\n");

	// Second line: Stem1-R, IS, Stem2-R
	//printf ("                    ");
	fprintf (output_file, "                    ");
	for (i=0; i<knot.s1; i++) {
	 //printf("%c", list[A5-i]);
	  switch (pcase) {
	    case 0:
	      fprintf (output_file, "%c", toupper(list[A5-i]));
	      break;
	    case 1:
	      fprintf (output_file, "%c", toupper(list[A5-i+bulge]));
	      break;
	    case 2:
	      fprintf (output_file, "%c", toupper(list[A5-i+bulge]));
	      break;
	    case 3:
	      if ((A5-i)>knot.blg)
	        fprintf (output_file, "%c", toupper(list[A5-i+bulge]));
	      else if ((A5-i)== knot.blg) {
	        fprintf (output_file, "*");
	        fprintf (output_file, "%c", toupper(list[A5-i+bulge]));
		fprintf (output_file, "%c", toupper(list[A5-i]));
	      } else 
		fprintf (output_file, "%c", toupper(list[A5-i]));
	      break;
	    case 4:
	      fprintf (output_file, "%c", toupper(list[A5-i]));
	      break;
  
	  } // switch
	}

       //printf("  ");
        fprintf (output_file, "  ");

	for (i=1; i<=knot.is; i++) {
	 //printf("%c", list[A4-i]);
	  switch (pcase) {
	    case 0:
	    case 3:
	    case 4:
	      fprintf (output_file, "%c", toupper(list[A4-i]));
	      break;
	    case 1:
	    case 2:
	      fprintf (output_file, "%c", toupper(list[A4-i+bulge]));
	      break;
	  } // switch
	}

       //printf("  ");
        fprintf (output_file, "  ");

	for (i=0; i<knot.s2; i++) {
	 //printf("%c", list[A3-i]);
	  switch (pcase) {
	    case 0:
	    case 3:
	    case 4:
	      fprintf (output_file, "%c", toupper(list[A3-i]));
	      break;
	    case 1:
	      fprintf (output_file, "%c", toupper(list[A3-i+bulge]));
	      break;
	    case 2:
	      if ((A3-i)>knot.blg)
	        fprintf (output_file, "%c", toupper(list[A3-i+bulge]));
	      else if ((A3-i)== knot.blg) {
	        fprintf (output_file, "*");
	        fprintf (output_file, "%c", toupper(list[A3-i+bulge]));
		fprintf (output_file, "%c", toupper(list[A3-i]));
	      } else 
		fprintf (output_file, "%c", toupper(list[A3-i]));
	      break;
	  } // switch
	}

       //printf("\n\n");
        fprintf (output_file, "\n\n");

	// Third line: Pre-Stem1, Stem1-L, Stem2-L
	for (i=20; i>0; i--) {
	  if ((A0-i)>=0) {
	   //printf("%c", list[A0-i]);
	    fprintf (output_file, "%c", toupper(list[A0-i]));
	  } else {
	   //printf(" ");
	    fprintf (output_file, " ");
	  }
	}

	for (i=0; i<knot.s1; i++) {
	 //printf("%c", list[A0+i]);
	  switch (pcase) {
	    case 0:
	    case 2:
	    case 3:
	    case 4:
	      fprintf (output_file, "%c", toupper(list[A0+i]));
	      break;
	    case 1:
	      if ((A0+i)>knot.blg)
	        fprintf (output_file, "%c", toupper(list[A0+i+bulge]));
	      else if ((A0+i)== knot.blg) {
	        fprintf (output_file, "*");
	        fprintf (output_file, "%c", toupper(list[A0+i]));
		fprintf (output_file, "%c", toupper(list[A0+i+bulge]));
	      } else 
		fprintf (output_file, "%c", toupper(list[A0+i]));
	      break;
	  } // switch
	}

       //printf("  ");
        fprintf (output_file, "  ");

	for (i=1; i<=knot.is; i++) {
	 //printf(" ");
	  fprintf (output_file, " ");
	}

       //printf("  ");
        fprintf (output_file, "  ");
	
	for (i=0; i<knot.s2; i++) {
	  switch (pcase) {
	    case 0:
	      fprintf (output_file, "%c", toupper(list[A6+i]));
	      break;
	    case 1:
	    case 2:
	    case 3:
	      fprintf (output_file, "%c", toupper(list[A6+i+bulge]));
	      break;
	    case 4:
	      //printf("%c", toupper(list[A6+i]);
	      if ((A6+i)>knot.blg)
	        fprintf (output_file, "%c", toupper(list[A6+i+bulge]));
	      else if ((A6+i)== knot.blg) {
	        fprintf (output_file, "*");
	        fprintf (output_file, "%c", toupper(list[A6+i]));
		fprintf (output_file, "%c", toupper(list[A6+i+bulge]));
	      } else 
		fprintf (output_file, "%c", toupper(list[A6+i]));
	      break;
	  } // switch
	}

	// 2 more after A7
	for (i=1; i<3; i++) {
	  fprintf (output_file, "%c", toupper(list[A7+i]));
	}

       //printf("\n\n");	
        fprintf (output_file, "\n\n");	
}

void phyloRna_c::timer(int sec, int msec) {
  fprintf (output_file, "Time taken to run: %d seconds %d milliseconds\n", sec, msec);
  fprintf(output_file, "-------------------------------------------------\n");
}

char *phyloRna_c::get_next_line (FILE *fp) {
	// try to allocate some memory
	char *buf = (char *)malloc (INIT_BUFF_SIZE);
	
	// return NULL on failure
	if (buf == NULL) {
		return NULL; }
	
	char c;
	int bytes_read = 0;
	
	// get the line from the file
	do {
		// get the next char in the file
		c = getc (fp);
		bytes_read++;

		// EOF or line feed or carriage return gets out of the loop
		if (c == EOF || c == '\n' || c == '\r') {
			continue; }

		buf[bytes_read - 1] = c;
		
		// that's all we have room for
		// output a warning and return the line so far
		if (bytes_read == (INIT_BUFF_SIZE - 2)) {
			printf ("-------------------------------------------------\n");
			printf ("get_next_line_stdio: line exceeded INIT_BUFF_SIZE\n");
			printf ("-------------------------------------------------\n");
			
			// get out of here
			break; }
	} while (c != EOF && c != '\n' && c != '\r');

	// major fail
	if (c == EOF && bytes_read == 1) {
		// have to free this since it was allocated
		free (buf);

		return NULL; }

	// append a terminator to the end of the new string
	buf[bytes_read - 1] = '\0';

	// return the line buffer
	return buf;
}


#endif //_phylo_h
