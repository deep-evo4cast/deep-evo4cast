All commands were executed on a Debian server with 2 Intel(R) Xeon(R)
Gold 6248R CPU @ 3.00GHz CPUs offering 48 cores and 96 threads, with
1TB of RAM using a fish shell and a python3.9.2 interpreter with
biopython and pyre2 packages installed. The amount of RAM may be
crucial for some steps below.

The spikeprot0125.tar file downloaded from GISAID was filtered to
remove all partial sequences, with length less than 1260, using the
filter-1260.py script, generating a "spike-1260.fasta2" FASTA file (in
FASTA2 BioPython format).

	python3 filter-1260.py

The last 18 amino acids of our predicted RBDs were extracted using
Unix filters cut, sort and uniq as follows:

	cut -b10-27 ddg_mono385_betterdg_noCYS.txt | sort | uniq > last18

To generate search patterns for Finite State Automata based search in
GISAID, we inserted single wildcards characters '.' (allowing for one
match) for each gap in our 27 AAs binding motifs. This was done using
the "generate-patterns.py" script, producing the 'last18.pats' file.

	python3 generate-patterns.py 18 last18

This 'last18.pats' file was cut into 50 pieces using:

	split -n l/50 last18.pats -d last18-pats.

producing files 'last18-pats.??' with 2 digits from 00 to 49 at the
end. These files were analyzed in parallel using the following fish
shell commands:

	for i in (seq -w 00 49)
   	   python3 find-disjunct.py $i &
	end

The script searches for the disjunction of all the patterns in the
filtered GISAID FASTA file, producing last18-pats.??.hits files. Only
files numbered 30, 36, 39 showed hits. To collect all patterns in
disjunctions with some GISAID hits, the following command was used:

	cat (find . -maxdepth 1 -name '*.hits' -type f ! -size 0 |sed s/\\.hits//)

The output of this command was redirected to a new last18-pats file in
a new directory named iter2 and the same process was applied (a
symbolic link to the spike-1260.fasta2 was made in the new directory)
and the same splitting and search commands used:

	split -n l/50 last18.pats -d last18-pats.
	for i in (seq -w 00 49)
           python3 ../find-disjunct.py $i &
        end

The same process was applied with an "iter3" directory: a total of 102
motifs with hits was idetified. For the last iteration, this file was
splitted in 102 files with one motif per file and the same process
applied:

        split -n l/102 last18.pats -d last18-pats.
        for i in (seq -w 000 101)
           python3 ../find-disjunct.py $i &
        end

This gave a final list of 9 patterns with hits that were collected in
a file named "last18-match". The gap characters were removed using:

	sed -e s/\\.//g last18-match > last18-match-dense 

and the initial list of putative RBMs filtered for occurrence of the
motifs:

	grep -f last18-match-dense ../ddg_mono385_betterdg_noCYS.txt | cut -b1-27 > last18-matches

This resulted in a list of 102 putative RBMs with potential matches in
GISAID.

To generate search regexps, we inserted single wildcards characters
'.' (allowing for one match) for each gap in our 27 AAs binding
motifs. This was done using the "generate-patterns.py" script,
producing the 'last18-matches.pats' file.

	python3	 generate-patterns.py 27 last18-matches

The presence of any of these motifs in the GISAID file was tested using:

	grep -f last18-matches.pats spike-1260.fasta2

showing no hit. 



We therefore decided to look for hits at Hamming distance 1. These
would be either putative RBMs with a mismatch in their first 9 amino
acids and a perfect match on the last 18 amino acids (category A) or
putative RBMs with a mismatch in the last 18 amino acids (category B).



For category A, we used the set of 102 RBMs in last18-matches.pats
(that have a perfect match on their last 18 characters) and inserted a
'.' in each non gap position, from position 1 to 9 in the initial
ungapped RBM. The following command:

	python3 ../generate-wc-sed.py
	
generated a list of 27 sed commands inserting a wildcard character at
each successive non-gap position in the motif. Only the 9 first
commands were executed:

	mkdir A-404 A-405 A-417 A-445 A-447 A-449 A-453 A-455 A-456
	sed s/././1 last18-matches.pats > A-404/last18-matches-wc-404.pats
   sed s/././2 last18-matches.pats > A-405/last18-matches-wc-405.pats
   sed s/././14 last18-matches.pats > A-417/last18-matches-wc-417.pats
   sed s/././42 last18-matches.pats > A-445/last18-matches-wc-445.pats
   sed s/././44 last18-matches.pats > A-447/last18-matches-wc-447.pats
   sed s/././46 last18-matches.pats > A-449/last18-matches-wc-449.pats
   sed s/././50 last18-matches.pats > A-453/last18-matches-wc-453.pats
   sed s/././52 last18-matches.pats > A-455/last18-matches-wc-455.pats
   sed s/././53 last18-matches.pats > A-456/last18-matches-wc-456.pats

In each of the A-404 A-405 A-417 A-445 A-447 A-449 A-453 A-455 A-456
folders, the find-disjunct2.py script, that uses re2 as before but
with a larger search window was used to check for matches using:

	for d in A-*
	   cd $d; ln -s ../spike-1260.fasta2
	   python3 ../find-disjunctA.py (string sub --start=-3 $d) &
	   cd ..
	end

showing that only mismatches at positions 405 and 445 could generate
hits. In these 2 cases, the 102 patterns were split in one line files
as above and the same search done. As an example, the commands
executed for position 405 were:

	split -n l/102 last18-matches-wc-405.pats -d last18-matches-wc-405-pats.
	for i in (seq -w 000 101)
            mkdir $i
            mv last18-matches-wc-405-pats.$i $i/last18-matches-wc-405.pats
	    cd $i; ln -s ../spike-1260.fasta2
	    python3 ../../find-disjunctA.py 405 &
	    cd ..
        end

This allowed to identify 1 putative RBM with a mismatch at position
405 (file 028) and 7 with a mismatch at position 445 (files 068, 091,
092, 093, 094, 095, 096) that had hits in GISAID (and the
corresponding matches).



For category B, we created a 18.pos file (provided) that contains the
positions of the last 18 amino acids in our RBM motif and created one
directory per position:

	for i in (cat 18.pos)
	   mkdir 18-$i
	   cd 18-$i
	   ln -s ../spike-1260.fasta2
	   cd ..
	end


From the last18.pats file above (containing 152487 motifs with gaps)
and for each non-gap position, we inserted a wildcard "." character at
each line, and removed duplicate lines using sed commands generated by
the generate-wc-sed18.py script:

	python3 generate-wc-sed18.py

which produced

	sed s/././1 last18.pats |sort | uniq > 18-473/last18.pats
	sed s/././3 last18.pats |sort | uniq > 18-475/last18.pats
	sed s/././4 last18.pats |sort | uniq > 18-476/last18.pats
	sed s/././12 last18.pats |sort | uniq > 18-484/last18.pats
	sed s/././14 last18.pats |sort | uniq > 18-486/last18.pats
	sed s/././15 last18.pats |sort | uniq > 18-487/last18.pats
	sed s/././17 last18.pats |sort | uniq > 18-489/last18.pats
	sed s/././21 last18.pats |sort | uniq > 18-493/last18.pats
	sed s/././22 last18.pats |sort | uniq > 18-494/last18.pats
	sed s/././23 last18.pats |sort | uniq > 18-495/last18.pats
	sed s/././24 last18.pats |sort | uniq > 18-496/last18.pats
	sed s/././26 last18.pats |sort | uniq > 18-498/last18.pats
	sed s/././28 last18.pats |sort | uniq > 18-500/last18.pats
	sed s/././29 last18.pats |sort | uniq > 18-501/last18.pats
	sed s/././30 last18.pats |sort | uniq > 18-502/last18.pats
	sed s/././31 last18.pats |sort | uniq > 18-503/last18.pats
	sed s/././32 last18.pats |sort | uniq > 18-504/last18.pats
	sed s/././33 last18.pats |sort | uniq > 18-505/last18.pats

The rest of the search was performed as in the initial search (with
motifs containing one wildcard instead of exact motifs), until a set
of patterns with hits could be identified, for each wildcard.

For each position, a list of 8 to 43 patterns were identified. Using
the remove-gaps18.py, wildcards characters corresponding to gaps were
removed from each pattern, leaving lists of patterns with only 1
wildcard.

These patterns were used to scan the list of all putative
RBMs (ddg_mono385_betterdg_noCYS.txt) using grep. Lists of 102 to
25870 paterns were collected.

As previously, wildcard "." were inserted at gap positions, using the
generate-patterns.py script. For each position, an original wildcard
character matching the position was added (using the
generate-wc-sed-final.py script) and the existence of hits with any of
the pattern was tested using the find-disjunctB-full.py.

Only wildcard positions 494, 498, 503 and 504 showed GISAID hits. Each
case was further dissected using the same splitting strategy as above
(using the find-disjunctB-last.py).

