// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * pkscan.c
// *
// * last modified: Sep-30-2023
// * Modified to process only in 5' region
// *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

//#include "stdafx.h"
#include <iostream>
#include "pkscan.h"

int main (int argc, char *argv[]) {
	// check for valid command line arguments
	time_t end, start;
	float time;
        bool verbose = false;
        
	start = clock();
	if (argc < 2) {
		fprintf (stderr, "\nusage: phylo [file_name] ... \n");
		exit (EXIT_FAILURE); 
	}
	char* filename;
	filename = argv[1];

	int stem11, stem12, stem21, stem22, lp11, lp12, lp21, lp22, is11, is12;
	stem11 = S1_MIN;
	stem12 = S1_MAX;
	stem21 = S2_MIN;
	stem22 = S2_MAX;
	lp11 = L1_MIN;
	lp12 = L1_MAX;
	lp21 = L2_MIN;
	lp22 = L2_MAX;
	is11 = IS_MIN;
	is12 = IS_MAX;

	if (argc >= 2) {
          int i=2;
	  while (i < argc) { /* We will iterate over argv[] to get the parameters stored inside.
			      * Note that we're starting on 1 because we don't need to know the 
			      * path of the program, which is stored in argv[0] */
            if (i + 1 != argc) // Check that we haven't finished parsing already
	      printf("parameter (%d of %d) %s!\n", i, argc, argv[i]);
            
	    if (strcmp(argv[i], "-s1") == 0) {
                    // We know the next argument *should* be the filename:
                    stem11 = atoi(argv[i + 1]);
		    stem12 = atoi(argv[i + 2]);
		    i = i+3;
	    } else if (strcmp(argv[i], "-s2") == 0) {
                    stem21 = atoi(argv[i + 1]);
		    stem22 = atoi(argv[i + 2]);
		    i = i+3;
	    } else if (strcmp(argv[i],  "-l1") == 0) {
                    lp11 = atoi(argv[i + 1]);
		    lp12 = atoi(argv[i + 2]);
		    i = i+3;
	    } else if (strcmp(argv[i], "-l2") == 0) {
                    lp21 = atoi(argv[i + 1]);
		    lp22 = atoi(argv[i + 2]);
		    i = i+3;
	    } else if (strcmp(argv[i], "-is") == 0) {
                    is11 = atoi(argv[i + 1]);
		    is12 = atoi(argv[i + 2]);
		    i = i+3;
            } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {                    
              // Check if "-v" or "--verbose" option is provided as a command-line argument
              verbose = true;
              i=1+2;
	    } else {
		  printf("Invalid parameter (of %d) (%d, %s)!\n", argc, i, argv[i]);
		  exit(0);
            }
           
	  } // while
	} // if (argc>3)

	// initialize the protein database class using the
	// command line arguments to set the file_name and max_dist
	pkScan_c  seq(argv[1], stem11, stem12, stem21, stem22,
		lp11, lp12, lp21, lp22, is11, is12);

        // set verbose mode
        seq.setVerboseMode(verbose);
        
	// load the data we're interested in
	seq.load_data ();

	// Get rid of possible base-pairing
	seq.find_pseudoknot();

	//	seq.load_result();

	// 
	seq.print_engy_knots ();


	// 
	//seq.print_knots ();

	// write files with indirect nucleotide to hetatm to amino acid connections
	//pdb.connect_nuc_het_ami ();

	end = clock();
	time = (end-start)/float(CLOCKS_PER_SEC)*1000;
	int minute=0, sec=0, msec=0;
	if (time > 1000) {
		sec = time/1000;
		msec = time - sec*1000;
	}
	printf("Time taken to run: %d seconds %d milliseconds\n", sec, msec);

	seq.timer(sec, msec);

		//hit any key to end
	//	cout << "Hit any key to terminate this program\n";
	//cin.get();
	// smooth sailing
	return 0;	
}
