// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * pkscan.h
// *
// * last modified: Sep-30-2023
// * Modified to process only in 5' region
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
#define S1_MIN       5     //5, 3
#define S1_MAX       15    //20
#define S2_MIN       5     //5, 3
#define S2_MAX       10
#define L1_MIN       1
#define L1_MAX       10    //5
#define L2_MIN       5
#define L2_MAX       35    //35
#define IS_MIN       0
#define IS_MAX       0
#define MIN_MIS_LEN  6
#define MIN_SUM_LEN  8
#define INIT_BUFF_SIZE		200
#define ENGY_MAX       1000
#define CHECK_MAX    30
//#define SS_MIN       3     //stem-stem
//#define SS_MAX       15     //stem-stem
//#define SL_MIN       4     //stem-loop
//#define SL_MAX       8     //stem-loop


typedef struct {
        int type;
        int num;
        double engy;
        double delta;
	int start;
        int end;
	int s1;
	int s2;
	int l1;
	int l2;
	int is;
        int mismatch;
        int  mis_s1;
        int  mis_s2;
        int blg;
} pseudoknot_t;

// class definition * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class pkScan_c {
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
	double max_dist;
	vector <pseudoknot_t> knots;
	int total;
	int allrank;
	FILE *output_file;
	int bulge;
	int mismatch;
	pseudoknot_t knot_engy[ENGY_MAX];
  int stloop;
  bool verboseMode; // 0 - simplify; 1 - detail

public:
  static const int SS_MIN = 5;     //stem-stem (3)
  static const int SS_MAX = 10;     //stem-stem
  static const int SL_MIN = 3;     //stem-loop (4)
  static const int SL_MAX = 8;     //stem-loop
  static const float MIN_ENGY;
  static const int Frt_Gap = 1;
  static const int End_Gap =2;
  static const int SL_L2_TAIL = 3;     //tail of remaining long L2 as new L2
  //static int stloop;
    
	pkScan_c ();
	pkScan_c (char *filename, int stem11, int stem12,
						int stem21, int stem22,
						int lp11, int lp12, int lp21, int lp22,
						int is11, int is12);
	~pkScan_c ();

	char *get_next_line (FILE *fp);

  void stloop_inc();
  int get_stloop();
	int load_data ();
	int load_result();
	void find_pseudoknot();
	bool is_pseudoknot(pseudoknot_t &knot);
	bool base_pair (int p1, int p2);
	void print_knot (pseudoknot_t knot);
        void find_stem_loop(pseudoknot_t knot);
	void timer(int sec, int msec);
	pseudoknot_t extract_result_from_line (char *line);
	pseudoknot_t cal_energy(pseudoknot_t knot);
	bool is_GU_pair(int p1, int p2, int p3, int p4);
	double get_wc_engy(int p1, int p2);
	double get_gu_engy(int p1, int p2, int p3, int p4);
	double get_3d_ends_engy(int p1, int p2, int p3);
	double get_5d_ends_engy(int p1, int p2, int p3);
	int comp_store_engy(pseudoknot_t knot);
  void insert_sort_engylist(pseudoknot_t knot, int pk_num);
	bool comp_to_keep(pseudoknot_t knot);
	void print_engy_knots();
  void setVerboseMode(bool verbose);
  void print_apical_loop(int A5, int A6, int l2_len, int st_pair, int l2_tail);  
};

const float pkScan_c::MIN_ENGY =20.0;

// methods  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


// constructor
pkScan_c::pkScan_c () {
	// to be safe, we'll initialize a couple things
	file_name = NULL;
}


// constructor that allows us to set the vlaues for file_name and max_dist
pkScan_c::pkScan_c (char *filename, int stem11, int stem12,
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

	stloop = 0;
	output_file = fopen (fname, "wt");
	  if (output_file == NULL) {
	    printf ("could not open %s for writing\n", fname);
	  }
	bulge = 0;
	mismatch = 0; // percent
}


// destructor
// clear all the lines that were loaded from the pdb file
pkScan_c::~pkScan_c () {
	// free all the original line buffers for the amino acids
  fclose (output_file);
}

void pkScan_c::stloop_inc() {
  stloop++;
}

int pkScan_c::get_stloop() {
  return stloop;
}

void pkScan_c::setVerboseMode(bool verbose) {
  verboseMode = verbose;
}

// load all the sequences from the input file
// assumes all atoms in a sequence will be listed consecutively
int pkScan_c::load_data (void) {
  string str;
  string acc ="";
  int l;
  int i;
  
	// let the user know where we are and what we're doing
	printf ("Sequence #: %s \n", file_name);
	fprintf (output_file, "Sequence #: %s \n", file_name);
        fflush(output_file); // Flush the output buffer
        
	FILE *file_pointer;			// pointer the the input file

	// try to open the input file
	file_pointer = fopen (file_name, "rt");

	// bail if we failed to open the file
	if (file_pointer == NULL) {
		printf ("error opening file: %s\n", file_name);
		return -1; }

	char *line;  // the line loaded from the input file
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
	fprintf(output_file, "Start=%d End=%d\n", cds_start, cds_end);
        fflush(output_file); // Flush the output buffer
        
        //l=0;
	while (1) {
		// get the next line from the database file
		line = get_next_line (file_pointer);
                //l++;
		// check for EOF
		if (line == NULL) {
                   break;
                }
		//fprintf (output_file, "%d: %s\n", l, line);
                //fflush(output_file); // Flush the output buffer
                
		for (i = 0; line[i] != '\0'; i++) {
			//cout << str.at(j) << "*";
			list.push_back(line[i]);
		}
		// read the nt until we reach the line where cds_start is on
		//if (list.size() >= cds_start) break;

	}

	// close the input file
	fclose (file_pointer);

	// tell us about what just happened
        //	printf (" (%lu) nucleotides\n", list.size ());
        fprintf (output_file, "loaded %lu nucleotides, cds_start %d to cds_end %d\n", list.size (), cds_start, cds_end);
        fflush(output_file); // Flush the output buffer
        
	//fprintf (output_file, "%s has (%lu) nucleotides\n", file_name, list.size ());

	// no problems encountered

	return 0;
}


void pkScan_c::find_pseudoknot (void) {
	int s1, s2;
	int l1, l2;
	int is;
	int loc, loc_max;
	bool is_knot;
	pseudoknot_t knot;
	bool found;
	int seq_max;
        int found_stems;
        int rank_pk;
        
        found_stems = 0;
	total = 0;
	allrank = 0;
        rank_pk = 0;
	seq_max = list.size();
	loc_max = seq_max - 2*stem_11 - 2*stem_21 - lp_11 - lp_21 - is_11;
	loc = 0;
	while (loc <= loc_max) {
          found=false;
          // for (s1 = stem_12;  s1 >= stem_11; s1--) {
          //   if (found && (knot.mis_s1==0)) break;
          //   else found=false;
          //   if ((loc+2*s1+2*stem_21+lp_11+lp_21+is_11) > seq_max) continue;
          //   for (s2 = stem_22;  s2 >= stem_21; s2--) {
          //     if (found &&(knot.mis_s2==0)) break;
          //     else found=false;
          //     if ((loc+2*s1+2*s2+lp_11+lp_21+is_11) > seq_max) continue;				
          //     for (is = is_11;  is <= is_12; is++) {
          //       if (found) break;
          //       if ((loc+2*s1+2*s2+lp_11+lp_21+is) > seq_max) break;
          //       for (l1 = lp_11;  l1 <= lp_12; l1++) {
          //         if (found) break;
          //         if ((loc+2*s1+2*s2+l1+lp_21+is) > seq_max) break;
          //         for (l2 = lp_21;  l2 <= lp_22; l2++) {
          //           if (found) break;
          //           if ((loc+2*s1+2*s2+l1+l2+is) > seq_max) break;

          for (s1 = stem_12;  s1 >= stem_11; s1--) {
            found=false;
            found_stems = 0;
            if ((loc+2*s1+2*stem_21+lp_11+lp_21+is_11) > seq_max) continue;
            for (s2 = stem_22;  s2 >= stem_21; s2--) {
              if (found &&((s1+s2)) < found_stems) break;
              found=false;
              found_stems = 0;
              if ((loc+2*s1+2*s2+lp_11+lp_21+is_11) > seq_max) continue;				
              for (is = is_11;  is <= is_12; is++) {
                found=false;
                found_stems = 0;
                if ((loc+2*s1+2*s2+lp_11+lp_21+is) > seq_max) break;
                for (l1 = lp_11;  l1 <= lp_12; l1++) {
                  found=false;
                  found_stems = 0;
                  if ((loc+2*s1+2*s2+l1+lp_21+is) > seq_max) break;
                  for (l2 = lp_21;  l2 <= lp_22; l2++) {
                    found=false;
                    found_stems = 0;
                    if ((loc+2*s1+2*s2+l1+l2+is) > seq_max) break;
                    
                    // nt # counting start from 1, so loc and loc_end add 1
                    // fprintf(output_file, "Here: loc_end=%d loc=%d s1=%d s2=%d l1=%d l2=%d is=%d\n", 
                    //      loc+2*s1+2*s2+l1+l2+is,loc+1, s1, s2, l1, l2, is); 
                    // c++ nt # count start from 0 to get the sequence code, so start and end is 1 less
                    // than printing
                    knot.start = loc;
                    knot.s1 = s1;
                    knot.s2 = s2;
                    knot.l1 = l1;
                    knot.l2 = l2;
                    knot.is = is;
                    knot.end = loc+2*s1+l1+2*s2+l2+is-1;
                    knot.blg = 0;
                    knot.mis_s2 = 0;
                    
                    if (bulge != 0) {
                      //is_knot = is_pseudoknot_bulge(knot);
                      if (is_knot) {
                        found = true;
                        if (knot.mismatch == 0)
                          //loc = loc + s1 - stem_11 + 1;
                          loc = loc + 1;
                      }
                    } else {
                      // if (comp_to_keep(knot)) {
                        is_knot = is_pseudoknot(knot);
			
                        if (is_knot) { // store seq is from [0..xx]
                          knot = cal_energy(knot);
                          if (knot.engy >= MIN_ENGY) {
                            total ++;
                            allrank = rank_pk +1;
                            knot.num = total;
                            rank_pk = comp_store_engy(knot);
                            found_stems = knot.s1 + knot.s2;
                            knots.push_back(knot);
                            found = true;
                            if (verboseMode) {                            
                              fprintf(output_file, "Found %d/%d: loc=%d s1=%d s2=%d l1=%d l2=%d is=%d end=%d engy=%f\n",
                                      rank_pk,knot.num, loc+1, s1, s2, l1, l2, is,knot.end+1, knot.engy);
                            }
                          }

                          //printf("Found %d: loc=%d s1=%d s2=%d l1=%d l2=%d is=%d\n",
                          //	total, loc+1, s1, s2, l1, l2, is);
                          //fprintf(output_file, "Found %d: loc=%d s1=%d s2=%d l1=%d l2=%d is=%d end=%d mismatch=%0d(s1:%0d, s2:%0d)\n",
                          
                          //total, loc+1, s1, s2, l1, l2, is,knot.end+1, knot.mismatch, knot.mis_s1, knot.mis_s2);
 
                          //print_knot(knot);
                          //Calculate the energy
                          //if (comp_to_keep(knot)) {

                          //fprintf(output_file, "Found %d: loc=%d s1=%d s2=%d l1=%d l2=%d is=%d end=%d engy=%f\n",
                          
                          //total, loc+1, s1, s2, l1, l2, is,knot.end+1, knot.engy);
                          //fprintf(output_file, "rank=%d engy=%f \n", knot.num, knot.engy);

                          //}
                          // if (knot.mismatch == 0)
                            //loc = loc + s1 - stem_11 + 1;
                          // loc = loc + 1;
                          
                        }//if (is_knot)
                        //} //if (comp_to_keep(knot))
                    }
                  } // is
                } // l2
              } //l1
            } // s2
          } // s1
          loc++;
          //if (!found) loc ++;
	} // loc
}



/*
bool pkScan_c::base_pair (int p1, int p2) {

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
			else 
			  return false;
			break;		
	  case 'g':
	  case 'G':
			if (list[p2]=='c' || list[p2]=='C') 
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
*/

bool pkScan_c::base_pair (int p1, int p2) {

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


pseudoknot_t pkScan_c::cal_energy(pseudoknot_t knot) {
	int i;
	double energy;
	int A0, A1, A2, A3, A4, A5, A6, A7;
	
	A0 = knot.start;
	A1 = A0 + knot.s1-1;
	A2 = A1 + knot.l1+1;
	A3 = A2 + knot.s2-1;
	A4 = A3 + knot.is+1;
	A5 = A4 + knot.s1-1;
	A6 = A5 + knot.l2+1;
	A7 = A6 + knot.s2-1;

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

void pkScan_c::insert_sort_engylist(pseudoknot_t knot, int pk_num) {
  int i, j;
  pseudoknot_t tmp;
  int lp;

  //allrank - all found pks
  if ((pk_num+1) < ENGY_MAX) 
    lp = pk_num;
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
  /* for (j=0; j<lp; j++) { */
  /*  fprintf(output_file, "\n Comp_store_engy: current=%d knot.engy=%f  knot_engy[%d]=%f \n", i, knot.engy, j, knot_engy[j].engy); */
  /* } */
}


int pkScan_c::comp_store_engy(pseudoknot_t knot) {
  int i, j;
  pseudoknot_t tmp;
  int lp;
  int current_pk;

  //allrank - all found pks
  if (allrank < ENGY_MAX) 
    lp = allrank;
  else
    lp = ENGY_MAX;

  //knot_engy[lp-1].engy = 0;
  for (i=0; i<lp; i++) {
    if ((knot.start >= knot_engy[i].start) && (knot.end <= knot_engy[i].end)){
      if (knot.engy >= knot_engy[i].engy) {
        // remove knot_engy[i], insert then return
        for (j=i; j<(lp-1); j++) {
          knot_engy[j] = knot_engy[j+1];
        }
        insert_sort_engylist(knot,lp-1);
          if (allrank < ENGY_MAX) 
            current_pk = allrank - 1;
          else
            current_pk = ENGY_MAX;

        return current_pk;
      }
      // if less engy, continue with the i loop   
    } else if ((abs(knot.start - knot_engy[i].start) <= Frt_Gap)&&(abs(knot.end - knot_engy[i].end) <= End_Gap)){
      if (knot.engy > knot_engy[i].engy) {
        // remove knot_engy[i], insert then return
        for (j=i; j<(lp-1); j++) {
          knot_engy[j] = knot_engy[j+1];
        }
        insert_sort_engylist(knot,lp-1);
          if (allrank < ENGY_MAX) 
            current_pk = allrank - 1;
          else
            current_pk = ENGY_MAX;

        return current_pk;
      } // if
    } // if
  } // for
  // if get to this line haven't return, then we get an extra pk. Insert and sort this pk
  insert_sort_engylist(knot,lp);
  return lp;
  /* for (j=0; j<lp; j++) { */
  /*  fprintf(output_file, "\n Comp_store_engy: current=%d knot.engy=%f  knot_engy[%d]=%f \n", i, knot.engy, j, knot_engy[j].engy); */
  /* } */
}


bool pkScan_c::comp_to_keep(pseudoknot_t knot) {
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



void pkScan_c::print_engy_knots() {
  int i, j;
  int lp;
  int rank;
  bool keep;
  int A2, A5, A6, B2, B5, B6;
	

  //printf("knot_num = %d \n", knot_num);
  //fprintf(output_file,"\n\n-------------------- Energy Rank ------------------------------\n\n");
  if (allrank < ENGY_MAX) 
    lp = allrank;
  else
    lp = ENGY_MAX;

  rank = 1;
  
  if (verboseMode) {                            
    fprintf(output_file, "\n Collect %d pseudoknots! \n", lp);
    for (i=0; i<lp; i++) {
      fprintf(output_file, "%d/%d knot_energy: Energy = %f  start= %d end=%d \n", 
              i+1, knot_engy[i].num, knot_engy[i].engy, knot_engy[i].start+1,knot_engy[i].end+1);
    }
  }
  
  for (i=0; i<lp; i++) {
    /* fprintf(output_file, "\nEnergy %d: -%f: Start=%d S1=%d S2=%d L1=%d L2=%d L3=%d End=%d\n\n", */

    /* i, knot_engy[i].engy, knot_engy[i].start+1, knot_engy[i].s1, knot_engy[i].s2, knot_engy[i].l1, knot_engy[i].l2, knot_engy[i].is, knot_engy[i].end+1); */
    /* fprintf(output_file, "\n # %d Energy= -%f \n", i, knot_engy[i].engy); */
    keep = true;
    if (knot_engy[i].engy < MIN_ENGY) continue;
    for (j=0; j<i; j++) {
       if ((knot_engy[i].start >= knot_engy[j].start)&&(knot_engy[i].end <= knot_engy[j].end)&&(knot_engy[i].engy <= knot_engy[j].engy)) {
	 A2 = knot_engy[i].start + knot_engy[i].s1 +knot_engy[i].l1;
	 A5 = A2 + knot_engy[i].s2+knot_engy[i].is+knot_engy[i].s1-1;
	 A6 = A5 + knot_engy[i].l2+1;
	 B2 = knot_engy[j].start + knot_engy[j].s1 +knot_engy[j].l1;
	 B5 = B2 + knot_engy[j].s2+knot_engy[j].is+knot_engy[j].s1-1;
	 B6 = B5 + knot_engy[j].l2+1;
         // assuming engy rank from high to low; there possible that stem1 first base-pair split or
         // stem2 last base-pair split 
	 if ((B2<=A2)&&(A5<=B5)&&(A6>=B6)) {
           if (verboseMode) {                            
             fprintf(output_file, "\n Remove #%d: Engy=%f (R:%d,%d) (S:%d,%d) i=%d A2,A5,A6 j=%d (R:%d,%d) B2,B5,B6: %d %d %d, %d %d %d (S:%d,%d) \n", 
                     knot_engy[i].num, knot_engy[i].engy, knot_engy[i].start+1,knot_engy[i].end+1,knot_engy[i].s1,knot_engy[i].s2,i+1,
                     j+1, knot_engy[j].start+1,knot_engy[j].end+1, A2, A5,A6,B2, B5,B6,knot_engy[j].s1,knot_engy[j].s2);
           }
	   keep = false;
	   break;
	 }
       }
    }// j
    if (!keep) continue;
    //	fprintf(output_file, "Energy %d: -%f: org_rank=%d start=%d s1=%d s2=%d l1=%d l2=%d is=%d end=%d\n",
    //		i+1, knot_engy[i].engy, knot_engy[i].num, knot_engy[i].start+1, knot_engy[i].s1, knot_engy[i].s2, knot_engy[i].l1, knot_engy[i].l2, knot_engy[i].is, knot_engy[i].end+1);
    if (verboseMode) {                            
      fprintf(output_file, "\n Keep #%d: rank=%d Energy = %f  start= %d end=%d i=%d A2,A5,A6 : %d %d %d s1=%d s2=%d\n", 
              knot_engy[i].num,rank, knot_engy[i].engy, knot_engy[i].start+1,knot_engy[i].end+1,i+1, A2, A5,A6,knot_engy[i].s1,knot_engy[i].s2);
    }
    fprintf(output_file, "\nEnergy %d: -%f: Start=%d S1=%d S2=%d L1=%d L2=%d L3=%d End=%d\n\n",

    rank, knot_engy[i].engy, knot_engy[i].start+1, knot_engy[i].s1, knot_engy[i].s2, knot_engy[i].l1, knot_engy[i].l2, knot_engy[i].is, knot_engy[i].end+1);
    print_knot(knot_engy[i]);
    //fprintf(output_file, "\n Searching for stem loop. \n");
    find_stem_loop(knot_engy[i]);
    rank += 1;
  }
}



bool pkScan_c::is_GU_pair(int p1, int p2, int p3, int p4) { // p1-p2:p3-p4 (p3=p1+1, p4=p2-1)

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

double pkScan_c::get_wc_engy(int p1, int p2) { 
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

double pkScan_c::get_gu_engy(int p1, int p2, int p3, int p4) { // p1-p2:p3-p4 p2=p1+1, p4=p2-1
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


double pkScan_c::get_3d_ends_engy(int p1, int p2, int p3) {
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

double pkScan_c::get_5d_ends_engy(int p1, int p2, int p3) {
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


bool pkScan_c::is_pseudoknot (pseudoknot_t &knot) {
	int A0, A1, A2, A3, A4, A5, A6, A7;
	int i;

	if (knot.s1+knot.s2<MIN_SUM_LEN) return false;
	A0 = knot.start;
	A1 = A0 + knot.s1-1;
	A2 = A1 + knot.l1+1;
	A3 = A2 + knot.s2-1;
	A4 = A3 + knot.is+1;
	A5 = A4 + knot.s1-1;
	A6 = A5 + knot.l2+1;
	A7 = A6 + knot.s2-1;
	
	int mis_s1=0;
	int mis_s2=0;
	

	if ( !base_pair(A0, A5) ) return false;
	if ( !base_pair(A1, A4) ) return false;
	if ( !base_pair(A6, A3) ) return false;
	if ( !base_pair(A7, A2) ) return false;

	for (i=1; i < (knot.s1-1); i++) {
		if ( !base_pair(A0+i, A5-i) ) {
		  if ((mismatch!=0)&&(knot.s1>MIN_MIS_LEN)) { // only allow mismatch when stem >5
				mis_s1 ++;
				if (mis_s1 > mismatch) return false;
			} else { // no bulge or mismatch allowed
				return false;
			}
		}
	}

	for (i=1; i < (knot.s2-1); i++) {
		if ( !base_pair(A6+i, A3-i) ) {
		  if ((mismatch!=0)&&(knot.s2>MIN_MIS_LEN)) { // only allow mismatch when stem >5
				mis_s2 ++;
				if ((mis_s1+mis_s2) > mismatch) return false;
			} else { // no bulge or mismatch allowed
				return false;
			}
		}
	}

	knot.mismatch =  mis_s1+mis_s2;
	knot.mis_s1 = mis_s1;
	knot.mis_s2 = mis_s2;
	
	return true;
}

void pkScan_c::find_stem_loop(pseudoknot_t knot) {
  int i, mt;  //match location
  int A0, A1, A2, A3, A4, A5, A6, A7;
  int stem1,l2_len;
  int st_lp_min, st_lp_max;
  int l2tail_min, l2tail_max;
  int l2tail;
  int ss1, ss2; //5'stem, 3'stem
  int s3, loop3; //stem 3, loop 3
  int st_pair, st3_max; // actual stem pairs, max of stem3 pairs
  int cal_max;
  bool max_pairs_reach;
  
  A0 = knot.start;
  A1 = A0 + knot.s1-1;
  A2 = A1 + knot.l1+1;
  A3 = A2 + knot.s2-1;
  A4 = A3 + knot.is+1;
  A5 = A4 + knot.s1-1;
  A6 = A5 + knot.l2+1;
  A7 = A6 + knot.s2-1;
  
  stem1 = knot.s1;
  l2_len = knot.l2;
  st_lp_min = SS_MIN*2 + SL_MIN; //10

  // the loop2 length is not long enough to form stem loop
  if (l2_len < (st_lp_min+ SL_L2_TAIL)) { //14
    return;
  }  

  // calculate the L2 tail (the remain length from knot.l2 used to form the apical stem loop)
  l2tail_max = l2_len -st_lp_min;
  st_lp_max = SS_MAX*2 + SL_MAX; //28
  if (l2_len > (st_lp_max + SL_L2_TAIL)) //32
    l2tail_min = l2_len -st_lp_max;
  else
    l2tail_min = SL_L2_TAIL;
  
  //fprintf(output_file, "l2_len=%d l2tail_min=%d l2tail_max=%d \n", l2_len,l2tail_min, l2tail_max);
  for (l2tail=l2tail_min; l2tail<=l2tail_max; l2tail++){
    st_pair = 0;
    max_pairs_reach = false;
    
    ss1 = A5 + 1;
    ss2 = A6 - l2tail -1;

    if (!base_pair(ss1,ss2)) {
      continue;  // contine to next l2tail iteration
    }
    //Got first base-pair
    st_pair += 1; //
    cal_max = (l2_len - l2tail - SL_MIN)/2; // max possible pairs,the fractional part discarded
    st3_max = (cal_max > SS_MAX) ? (SS_MAX-1): (cal_max-1);  // -1 because we got one base-pair before coming into here
    for (s3=st3_max; s3>0; s3--){ // s3 used as iteration variable
      ss1 ++;
      ss2 --;
      if (!base_pair(ss1,ss2)) {
        if (st_pair < SS_MIN) {
          //fprintf(output_file, "Here l2_tail=%d s3=%d st_pair=%d !\n", l2tail, s3, st_pair);
          break; //break out of s3 loop and go to l2tail loop
        } else {
          // check if the stem loop (l3) satisfy min, max requirement
          loop3 = l2_len - l2tail - (st_pair*2);
          if ((loop3>=SL_MIN)&&(loop3<=SL_MAX)){
            // enough pairs, and loop3 is in range
            // found one and print stem loop
            // st_pair,(ABC)loop3(CBA)-l2tail
            max_pairs_reach = true;
            stloop_inc();
            print_apical_loop(A5, A6, l2_len, st_pair, l2tail);
          } // if loop3
          break; 
        } // if st_pair
      } else {
        st_pair += 1;
        //fprintf(output_file, "Here l2_tail=%d s3=%d st_pair=%d !\n", l2tail, s3, st_pair);
      } // if base_pair
    } // for s3
    // here capture the last one, which has the longest (SS_MAX) pairs,
    //and loop3 still in range (loop3>SL_MIN+1). We only need to check if loop3 exceed SL_MAX
    if ((!max_pairs_reach) && (st_pair >= SS_MIN)) { //never print out before
      loop3 = l2_len - l2tail - (st_pair*2);
      if ((loop3>=SL_MIN)&&(loop3<=SL_MAX)){
        // enough pairs, and loop3 is in range
        // found one and print stem loop
        // st_pair,(ABC)loop3(CBA)-l2tail
        max_pairs_reach = true;
        stloop_inc();
        print_apical_loop(A5, A6, l2_len, st_pair, l2tail);
      } // if loop3
    }
  } // for l2tail
    
}

void pkScan_c::print_apical_loop(int A5, int A6, int l2_len, int st_pair, int l2_tail) {
  int l3;
  int i;

  l3 =  l2_len - l2_tail - (st_pair*2);
  // st_pair,(ABC)AUGCAUGC(CBA)-REMAININGNUCLETIDE
  // st_pair,(ABC)l3(CBA)-l2_tail

  fprintf (output_file, "Found #%d stem loop within loop2: S3=%d AL=%d L2_tail=%d loop2(%d-%d) cds_start=%d\n",
           get_stloop(), st_pair, l3, l2_tail, A5+2, A6-1, cds_start);
  //for (i=A5+1;i<=A5+1+mt; i++){
  for (i=A5+1;i<=A6-1; i++){
    if ((i == (A5+1)) or (i == (A5+st_pair+l3+1))){
      fprintf (output_file, "(");
    }
    fprintf (output_file, "%c", toupper(list[i]));
    if ((i == (A5+st_pair)) or (i== (A6-l2_tail-1))) {
      fprintf (output_file, ")");
    }
  }
  fprintf (output_file, "\n\n");
  
}


// void pkScan_c::find_stem_loop(pseudoknot_t knot) {
//   int i, mt;  //match location
//   int stem_min, stem_max, ss1, ss2;
//   int stem1,l2_len;
//   int sp, st_pair;
//   int p1, p2;
//   int A0, A1, A2, A3, A4, A5, A6, A7;
//   int st_lp_min;
//   int st_1mt_min, st_1mt_max; //first match locates in the sequence
//   int tmp1;

//   A0 = knot.start;
//   A1 = A0 + knot.s1-1;
//   A2 = A1 + knot.l1+1;
//   A3 = A2 + knot.s2-1;
//   A4 = A3 + knot.is+1;
//   A5 = A4 + knot.s1-1;
//   A6 = A5 + knot.l2+1;
//   A7 = A6 + knot.s2-1;
  
//   stem1 = knot.s1;
//   l2_len = knot.l2;
//   st_lp_min = SS_MIN*2 + SL_MIN + stem1;
//   //st_pair = 0;
//   // the loop2 length is not long enough to form stem loop
//   if (l2_len < st_lp_min) {
//     //fprintf(output_file, "\n Loop2 is too short: %d, need to be at least %d \n", l2_len, st_lp_min);
//     return;
//   }  

//   //calculate the range of first nt of the stem loop
//   st_1mt_min = (SS_MIN - 1)+ SL_MIN + SS_MIN;  //2+4+3=9
//   st_1mt_max = (SS_MAX - 1)+ SL_MAX + SS_MAX;  //14+8+15=37
//   // check max of first match, to keep the # of nt long enough for covering the cpk stem1
//   if ((l2_len - st_1mt_max) < stem1) { st_1mt_max = l2_len - stem1; }
//   // The first of stem loop is A5+1
//   for (mt=st_1mt_min; mt<=st_1mt_max; mt++) {
//     ss1 = A5 + 1;
//     ss2 = ss1 + mt;
//     //fprintf(output_file, "A5=%d l2_len=%d ss1=%d ss2=%d mt=%d. \n", A5, l2_len,ss1, ss2, mt);
//     if (!base_pair(ss1,ss2)) {
//       //fprintf(output_file, "Base Pair: %c (%d),  %c (%d) does not match. \n", list[ss1], ss1, list[ss2], ss2);
//       continue;
//     }
//     //Got the match for the first nt, now from A5+1 to mt, we form a stem loop
//     //Check what is the max stem pairs we can form, with the ST_MIN
//     // C++ divison truncates the decimal part
//     st_pair =1;
//     stem_max = (ss2 - A5 - SL_MIN) / 2; // max possible pairs,the fractional part discarded
//     tmp1 = (ss2 - A5 - SL_MAX + 1) / 2; // round up, y/n => (y+n-1)/n 
//     stem_min = (tmp1>SS_MIN)?tmp1:SS_MIN; // min pairs to ensure AL loop length
//     //fprintf (output_file, "stem_max=%d stem_min=%d \n", stem_max, stem_min);
//     // the nucleotides is from A5+1 to A5+1+mt
//     while (st_pair < stem_max) {
//       ss1 ++;
//       ss2 --;
//       if (!base_pair(ss1,ss2)) {
// 	if (st_pair >= stem_min) {
// 	  // found one and print stem loop
// 	  // st_pair,(ABC)AUGCAUGC(CBA)-REMAININGNUCLETIDE
// 	  stloop_inc();
// 	  fprintf (output_file, "Found #%d stem loop within loop2: Start=%d, S3=%d AL=%d L2_tail=%d L2-range(0-%d)\n",
// 		   get_stloop(), A0+1, st_pair, mt-2*st_pair+1, l2_len-mt-1, mt);
// 	  //for (i=A5+1;i<=A5+1+mt; i++){
// 	  for (i=A5+1;i<=A6-1; i++){
// 	    if ((i == (A5+1)) or (i == (A5+1+mt-st_pair+1))){
// 		fprintf (output_file, "(");
// 	    }
// 	    fprintf (output_file, "%c", toupper(list[i]));
//  	    if ((i == (A5+st_pair)) or (i== (A5+1+mt))) {
// 		fprintf (output_file, ")");
// 	    }
// 	  }
// 	  fprintf (output_file, "\n\n");
// 	  //st_pair=0;
// 	} else {
// 	  //for debug purpose
// 	  //fprintf (output_file, "NOT Found! a stem loop of %d pairs:\n", st_pair);	  
// 	}
// 	break;
//       } else {  // if (!base_pair(ss1,ss2))
// 	st_pair ++;
//       }
//     } // while
//   } //for mt
  
// }


void pkScan_c::print_knot (pseudoknot_t knot) {
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
          if ((A7+i)<list.size())
            fprintf (output_file, "%c", toupper(list[A7+i]));
	}

       //printf("\n\n");	
        fprintf (output_file, "\n\n");	
}

void pkScan_c::timer(int sec, int msec) {
  fprintf (output_file, "Time taken to run: %d seconds %d milliseconds\n", sec, msec);
  fprintf(output_file, "-------------------------------------------------\n");
}

char *pkScan_c::get_next_line (FILE *fp) {
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
