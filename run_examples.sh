echo "./Achaz_stats -f 1 10 9 2 2 1 0 1 0 0 1 > ./fasta_checknsam10unfold_tests00.txt"
echo "..."
./Achaz_stats -f 1 10 9 2 2 1 0 1 0 0 1 > ./fasta_checknsam10unfold_tests00.txt 
echo "done"
echo "./Achaz_stats -f 0 10 10 2 2 2 0 > ./fasta_checknsam10fold_tests00.txt"
echo "..."
./Achaz_stats -f 0 10 10 2 2 2 0 > ./fasta_checknsam10fold_tests00.txt
echo "done"
echo "./Achaz_stats -f 1 10 50 0 10 5 4 3 2 1 0 > ./fasta_checknsam10unfold_tests01_full.txt"
echo "..."
./Achaz_stats -f 1 10 50 0 10 5 4 3 2 1 0 > ./fasta_checknsam10unfold_tests01_full.txt
echo "done"
echo "./Achaz_stats -s 1 10 1:50 3:10 4:5 5:4 6:3 7:2 8:1 > ./fasta_checknsam10unfold_tests01_part.txt"
echo "..."
./Achaz_stats -s 1 10 1:50 3:10 4:5 5:4 6:3 7:2 8:1 > ./fasta_checknsam10unfold_tests01_part.txt
echo "done"
echo "./Achaz_stats -d 25484 -s 1 5000 1:50 3:10 4:5 5:4 6:3 7:2 543:1 > ./fasta_checknsam10unfold_tests02_n5000.txt"
echo "..."
./Achaz_stats -d 25484 -s 1 5000 1:50 3:10 4:5 5:4 6:3 7:2 543:1 > ./fasta_checknsam10unfold_tests02_n5000.txt
echo "done"
echo "./Achaz_stats -d 25484 -s 1 4000 1:50 3:10 4:5 5:4 6:3 7:2 543:1 > ./fasta_checknsam10unfold_tests02_n4000.txt"
echo "..."
./Achaz_stats -d 25484 -s 1 4000 1:50 3:10 4:5 5:4 6:3 7:2 543:1 > ./fasta_checknsam10unfold_tests02_n4000.txt
echo "done"
echo "./Achaz_stats -d 25484 -s 1 4001 1:50 3:10 4:5 5:4 6:3 7:2 543:1 > ./fasta_checknsam10unfold_tests02_n4001.txt"
echo "..."
./Achaz_stats -d 25484 -s 1 4001 1:50 3:10 4:5 5:4 6:3 7:2 543:1 > ./fasta_checknsam10unfold_tests02_n4001.txt
echo "done"
echo "./Achaz_stats -d 46459 -s 1 10000 1:50 3:10 4:5 5:4 6:3 7:2 543:1 > ./fasta_checknsam10unfold_tests02_n10000.txt"
echo "..."
./Achaz_stats -d 46459 -s 1 10000 1:50 3:10 4:5 5:4 6:3 7:2 543:1 > ./fasta_checknsam10unfold_tests02_n10000.txt
echo "done"
echo "./Achaz_stats -d 88064 -s 1 1000000 1:50 3:10 4:5 5:4 6:3 7:2 > ./fasta_checknsam10unfold_tests02_n1e6.txt"
echo "..."
./Achaz_stats -d 88064 -s 1 1000000 1:50 3:10 4:5 5:4 6:3 7:2 > ./fasta_checknsam10unfold_tests02_n1e6.txt
echo "done"
echo "./Achaz_stats -d 88064 -s 1 1000000 1:50 3:10 4:5 5:4 6:3 7:2 543:1 > ./fasta_checknsam10unfold_tests02_n1e6b.txt"
echo "..."
./Achaz_stats -d 88064 -s 1 1000000 1:50 3:10 4:5 5:4 6:3 7:2 543:1 > ./fasta_checknsam10unfold_tests02_n1e6b.txt
echo "done"
echo "./Achaz_stats -f 1 20 77 38 8 0 0 9 0 9 8 0 8 0 0 0 0 0 0 0 0 > ./fasta_checknsam20unfold_test03.txt"
echo "..."
./Achaz_stats -f 1 20 77 38 8 0 0 9 0 9 8 0 8 0 0 0 0 0 0 0 0 > ./fasta_checknsam20unfold_test03.txt
echo "done"
echo "./Achaz_stats -f 0 3 1 > ./fasta_checknsam3fold_test04.txt"
echo "..."
./Achaz_stats -f 0 3 1 > ./fasta_checknsam3fold_test04.txt
echo "done"
echo "./Achaz_stats -f 0 4 1 0 > ./fasta_checknsam4fold_test04.txt"
echo "..."
./Achaz_stats -f 0 4 1 0 > ./fasta_checknsam4fold_test04.txt
echo "done"
echo "./Achaz_stats -f 0 5 1 0 > ./fasta_checknsam5fold_test05.txt"
echo "..."
./Achaz_stats -f 0 5 1 0 > ./fasta_checknsam5fold_test05.txt
echo "done"
echo "./Achaz_stats -f 0 5 0 1 > ./fasta_checknsam5fold_test06.txt"
echo "..."
./Achaz_stats -f 0 5 0 1 > ./fasta_checknsam5fold_test06.txt
echo "done"
echo "End Examples"
