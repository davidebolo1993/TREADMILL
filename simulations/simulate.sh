#!/bin/bash


if [ ! -f hs37d5.fa ]; then

	wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
	gzip -d hs37d5.fa.gz
	genome=$(readlink -f hs37d5.fa)

fi

error="85 90 95" #3 accuracy levels
coverage="20 50 100 200" #4 coverage levels

mkdir -p CGG_affected
cd CGG_affected

echo "Simulating affected samples"

#affected (expansion from 20 to 250 CGG motifs)

echo -e "X\t146993561\t146993628\ttandem repeat expansion\tCGG:230\t0" > hack.bed
echo -e "X\t146993561\t146993628" > treadmill.bed
chromosome=$(cut -f 1 hack.bed)
start=$(cut -f 2 hack.bed)
end=$(cut -f 3 hack.bed)
newstart=$((${start}-4000))
newend=$((${end}+4000))
echo -e ${chromosome}"\t"${newstart}"\t"${newend}"\t100.0\t100.0" > laser.bed
affinity="70.0" #this is overall fine if we want to cluster reads that belong to groups having large differences


if [ ! -f "${chromosome}.fa" ]; then

	samtools faidx ${genome} ${chromosome} > ${chromosome}.fa

fi

VISOR HACk -g ${chromosome}.fa -b hack.bed -o hack_out && cp ${chromosome}.fa hack_out/h2.fa

for j in ${error}; do

	mkdir -p ${j}_percent
	cd ${j}_percent
	echo "Simulating using identity " ${j}

	for l in ${coverage}; do

		echo "Simulating using coverage " ${l}
		mkdir -p ${l}X
		cd ${l}X

		if [ ! -f results_DBSCAN.txt ]; then

			echo -e "simulation\texpected\tpredicted\tcoverage\tidentity\talgorithm" > results_DBSCAN.txt
			echo -e "simulation\texpected\tpredicted\tcoverage\tidentity\talgorithm" > results_AGGLOMERATIVE.txt

			for i in {1..200}; do

				echo "Simulation " ${i}

				VISOR LASeR -g ../../${chromosome}.fa -s ../../hack_out -b ../../laser.bed -o ${i}_laser_out --coverage ${l} --length_mean 8000 --length_stdev 500 --tag --identity_min ${j}

				if [[ "${l}" == "20" ]]; then

					support="5"

				elif [[ "${l}" == "50" ]]; then

					support="10"

				elif [[ "${l}" == "100" ]]; then

					support="20"

				else

					support="40"

				fi

				TREADMILL RACE -fa ${genome} -bam ${i}_laser_out/sim.srt.bam -bed ../../treadmill.bed --motif CGG -o ${i}_treadmill_out_DBSCAN/${i}.bin --store --plot --support ${support} --similarity 90.0 --threads 4 --affinity ${affinity}
				TREADMILL TRAP -i ${i}_treadmill_out_DBSCAN/${i}.bin -o ${i}_treadmill_out_DBSCAN

				predicted=$(bcftools query -f '%RALN,%AL1N' ${i}_treadmill_out_DBSCAN/TREADMILL.vcf.gz)
				rn=$(bcftools query -f '%RN\n' ${i}_treadmill_out_DBSCAN/TREADMILL.vcf.gz)
				en=$((${rn} + 230))
				expected="$rn,$en"
				echo -e ${i}"\t"${expected}"\t"${predicted}"\t"${l}X"\t"${j}%"\tDBSCAN" >> results_DBSCAN.txt

				TREADMILL RACE -fa ${genome} -bam ${i}_laser_out/sim.srt.bam -bed ../../treadmill.bed --motif CGG -o ${i}_treadmill_out_DENDOGRAM/${i}.bin --similarity 90.0 --threads 4 --hierarchical_clustering --dendogram --store
				python /home/davide/tools/TREADMILL/scripts/clusteranalysis.py -d ${i}_treadmill_out_DENDOGRAM/${i}.bin -s ${i}_treadmill_out_DENDOGRAM/simmatrix.bin -o ${i}_treadmill_out_SILHOUETTE
				n_clusters=$(tail -n+2 ${i}_treadmill_out_SILHOUETTE/silhouettescores.tsv | cut -f 1 | head -1)
				TREADMILL RACE -fa ${genome} -bam ${i}_laser_out/sim.srt.bam -bed ../../treadmill.bed --motif CGG -o ${i}_treadmill_out_AGGLOMERATIVE/${i}.bin --similarity 90.0 --threads 4 --hierarchical_clustering --clusters ${n_clusters} --store
				TREADMILL TRAP -i ${i}_treadmill_out_AGGLOMERATIVE/${i}.bin -o ${i}_treadmill_out_AGGLOMERATIVE

				predicted=$(bcftools query -f '%RALN,%AL1N' ${i}_treadmill_out_AGGLOMERATIVE/TREADMILL.vcf.gz)
				rn=$(bcftools query -f '%RN\n' ${i}_treadmill_out_AGGLOMERATIVE/TREADMILL.vcf.gz)
				en=$((${rn} + 230))
				expected="$rn,$en"
				echo -e ${i}"\t"${expected}"\t"${predicted}"\t"${l}X"\t"${j}%"\tAGGLOMERATIVE" >> results_AGGLOMERATIVE.txt

			done

		else #file exists, only calculate TPR which is fast even if it does exist

			echo "Skipping simulations"

		fi
		
		#calculate TPR for the 2 methods

		echo -e "TPR\terrors\tcoverage\tidentity\talgorithm\tstatus" > TPR_DBSCAN.txt
		echo -e "TPR\terrors\tcoverage\tidentity\talgorithm\tstatus" > TPR_AGGLOMERATIVE.txt

		er="0 1 2 3 4 5 6 7 8 9 10" #allow up to 10 errors
		
		for z in ${er}; do

			TP=$(tail -n+2 results_DBSCAN.txt | awk -v var=${z} '{split($3,a,","); if (a[2] != "." && a[2]>=250-var && a[2]<=250+var) print "0"; else print "1"}' | grep -c "0") #0 are correct; 1 incorrect
			total=$(tail -n+2 results_DBSCAN.txt | grep -c "^")
			TPR=$(echo ${TP}/${total} | bc -l)
			echo -e ${TPR}"\t"${z}"\t"${l}X"\t"${j}%"\tDBSCAN\taffected" >> TPR_DBSCAN.txt

			TP=$(tail -n+2 results_AGGLOMERATIVE.txt | awk -v var=${z} '{split($3,a,","); if (a[2] != "." && a[2]>=250-var && a[2]<=250+var)  print "0"; else print "1"}' | grep -c "0") #0 are correct; 1 incorrect
			total=$(tail -n+2 results_AGGLOMERATIVE.txt | grep -c "^")
			TPR=$(echo ${TP}/${total} | bc -l)
			echo -e ${TPR}"\t"${z}"\t"${l}X"\t"${j}%"\tAGGLOMERATIVE\taffected" >> TPR_AGGLOMERATIVE.txt
		
		done
		cd ..
	done 
	cd ..
done
cd ..


genome="/home/davide/data/crispr_human/hs37d5.decoy.fa"
error="85 90 95" #3 error levels
coverage="20 50 100 200" #4 coverage levels

mkdir -p CGG_premutated
cd CGG_premutated

echo "Simulating premutated samples"

#premutated (expansion fro 20 to 150 CGG motifs)

echo -e "X\t146993561\t146993628\ttandem repeat expansion\tCGG:130\t0" > hack.bed
echo -e "X\t146993561\t146993628" > treadmill.bed
chromosome=$(cut -f 1 hack.bed)
start=$(cut -f 2 hack.bed)
end=$(cut -f 3 hack.bed)
newstart=$((${start}-4000))
newend=$((${end}+4000))
echo -e ${chromosome}"\t"${newstart}"\t"${newend}"\t100.0\t100.0" > laser.bed
affinity="70.0"


if [ ! -f "${chromosome}.fa" ]; then

	samtools faidx ${genome} ${chromosome} > ${chromosome}.fa

fi

VISOR HACk -g ${chromosome}.fa -b hack.bed -o hack_out && cp ${chromosome}.fa hack_out/h2.fa

for j in ${error}; do

	mkdir -p ${j}_percent
	cd ${j}_percent
	echo "Simulating using identity " ${j}

	for l in ${coverage}; do

		echo "Simulating using coverage " ${l}
		mkdir -p ${l}X
		cd ${l}X

		if [ ! -f results_DBSCAN.txt ]; then

			echo -e "simulation\texpected\tpredicted\tcoverage\tidentity\talgorithm" > results_DBSCAN.txt
			echo -e "simulation\texpected\tpredicted\tcoverage\tidentity\talgorithm" > results_AGGLOMERATIVE.txt

			for i in {1..200}; do

				echo "Simulation " ${i}

				VISOR LASeR -g ../../${chromosome}.fa -s ../../hack_out -b ../../laser.bed -o ${i}_laser_out --coverage ${l} --length_mean 8000 --length_stdev 500 --tag --identity_min ${j}

				if [[ "${l}" == "20" ]]; then

					support="5"

				elif [[ "${l}" == "50" ]]; then

					support="10"

				elif [[ "${l}" == "100" ]]; then

					support="20"

				else

					support="40"

				fi

				TREADMILL RACE -fa ${genome} -bam ${i}_laser_out/sim.srt.bam -bed ../../treadmill.bed --motif CGG -o ${i}_treadmill_out_DBSCAN/${i}.bin --store --plot --support ${support} --similarity 90.0 --threads 4 --affinity ${affinity}
				TREADMILL TRAP -i ${i}_treadmill_out_DBSCAN/${i}.bin -o ${i}_treadmill_out_DBSCAN

				predicted=$(bcftools query -f '%RALN,%AL1N' ${i}_treadmill_out_DBSCAN/TREADMILL.vcf.gz)
				rn=$(bcftools query -f '%RN\n' ${i}_treadmill_out_DBSCAN/TREADMILL.vcf.gz)
				en=$((${rn} + 130))
				expected="$rn,$en"
				echo -e ${i}"\t"${expected}"\t"${predicted}"\t"${l}X"\t"${j}%"\tDBSCAN" >> results_DBSCAN.txt

				TREADMILL RACE -fa ${genome} -bam ${i}_laser_out/sim.srt.bam -bed ../../treadmill.bed --motif CGG -o ${i}_treadmill_out_DENDOGRAM/${i}.bin --similarity 90.0 --threads 4 --hierarchical_clustering --dendogram --store
				python /home/davide/tools/TREADMILL/scripts/clusteranalysis.py -d ${i}_treadmill_out_DENDOGRAM/${i}.bin -s ${i}_treadmill_out_DENDOGRAM/simmatrix.bin -o ${i}_treadmill_out_SILHOUETTE
				n_clusters=$(tail -n+2 ${i}_treadmill_out_SILHOUETTE/silhouettescores.tsv | cut -f 1 | head -1)
				TREADMILL RACE -fa ${genome} -bam ${i}_laser_out/sim.srt.bam -bed ../../treadmill.bed --motif CGG -o ${i}_treadmill_out_AGGLOMERATIVE/${i}.bin --similarity 90.0 --threads 4 --hierarchical_clustering --clusters ${n_clusters} --store
				TREADMILL TRAP -i ${i}_treadmill_out_AGGLOMERATIVE/${i}.bin -o ${i}_treadmill_out_AGGLOMERATIVE

				predicted=$(bcftools query -f '%RALN,%AL1N' ${i}_treadmill_out_AGGLOMERATIVE/TREADMILL.vcf.gz)
				rn=$(bcftools query -f '%RN\n' ${i}_treadmill_out_AGGLOMERATIVE/TREADMILL.vcf.gz)
				en=$((${rn} + 130))
				expected="$rn,$en"
				echo -e ${i}"\t"${expected}"\t"${predicted}"\t"${l}X"\t"${j}%"\tAGGLOMERATIVE" >> results_AGGLOMERATIVE.txt

			done

		else 

			echo "Skipping simulations"

		fi 
		
		#calculate TPR for the 2 methods

		echo -e "TPR\terrors\tcoverage\tidentity\talgorithm\tstatus" > TPR_DBSCAN.txt
		echo -e "TPR\terrors\tcoverage\tidentity\talgorithm\tstatus" > TPR_AGGLOMERATIVE.txt

		er="0 1 2 3 4 5 6 7 8 9 10" #allow up to 10 errors

		for z in ${er}; do

			TP=$(tail -n+2 results_DBSCAN.txt | awk -v var=${z} '{split($3,a,","); if (a[2] != "." && a[2]>=150-var && a[2]<=150+var) print "0"; else print "1"}' | grep -c "0") #0 are correct; 1 incorrect
			total=$(tail -n+2 results_DBSCAN.txt | grep -c "^")
			TPR=$(echo ${TP}/${total} | bc -l)
			echo -e ${TPR}"\t"${z}"\t"${l}X"\t"${j}%"\tDBSCAN\tpremutated" >> TPR_DBSCAN.txt

			TP=$(tail -n+2 results_AGGLOMERATIVE.txt | awk -v var=${z} '{split($3,a,","); if (a[2] != "." && a[2]>=150-var && a[2]<=150+var)  print "0"; else print "1"}' | grep -c "0") #0 are correct; 1 incorrect
			total=$(tail -n+2 results_AGGLOMERATIVE.txt | grep -c "^")
			TPR=$(echo ${TP}/${total} | bc -l)
			echo -e ${TPR}"\t"${z}"\t"${l}X"\t"${j}%"\tAGGLOMERATIVE\tpremutated" >> TPR_AGGLOMERATIVE.txt
		
		done
		cd ..
	done 
	cd ..
done
cd ..

genome="/home/davide/data/crispr_human/hs37d5.decoy.fa"
error="85 90 95" #3 error levels
coverage="20 50 100 200" #4 coverage levels

mkdir -p CGG_greyzone
cd CGG_greyzone

echo "Simulating greyzone samples"

#greyzone (expansion from 20 to 50 CGG motifs)

echo -e "X\t146993561\t146993628\ttandem repeat expansion\tCGG:30\t0" > hack.bed
echo -e "X\t146993561\t146993628" > treadmill.bed
chromosome=$(cut -f 1 hack.bed)
start=$(cut -f 2 hack.bed)
end=$(cut -f 3 hack.bed)
newstart=$((${start}-4000))
newend=$((${end}+4000))
echo -e ${chromosome}"\t"${newstart}"\t"${newend}"\t100.0\t100.0" > laser.bed


if [ ! -f "${chromosome}.fa" ]; then

	samtools faidx ${genome} ${chromosome} > ${chromosome}.fa

fi

VISOR HACk -g ${chromosome}.fa -b hack.bed -o hack_out && cp ${chromosome}.fa hack_out/h2.fa

for j in ${error}; do

	mkdir -p ${j}_percent
	cd ${j}_percent
	echo "Simulating using identity " ${j}

	if [[ "${j}" == "85" ]]; then

		affinity="80"

	else

		affinity="85"
	
	fi

	for l in ${coverage}; do

		echo "Simulating using coverage " ${l}
		mkdir -p ${l}X
		cd ${l}X
		
		if [ ! -f results_DBSCAN.txt ]; then

			echo -e "simulation\texpected\tpredicted\tcoverage\tidentity\talgorithm" > results_DBSCAN.txt
			echo -e "simulation\texpected\tpredicted\tcoverage\tidentity\talgorithm" > results_AGGLOMERATIVE.txt

			for i in {1..200}; do

				echo "Simulation " ${i}

				VISOR LASeR -g ../../${chromosome}.fa -s ../../hack_out -b ../../laser.bed -o ${i}_laser_out --coverage ${l} --length_mean 8000 --length_stdev 500 --tag --identity_min ${j}

				if [[ "${l}" == "20" ]]; then

					support="5"

				elif [[ "${l}" == "50" ]]; then

					support="10"

				elif [[ "${l}" == "100" ]]; then

					support="20"

				else

					support="40"

				fi

				TREADMILL RACE -fa ${genome} -bam ${i}_laser_out/sim.srt.bam -bed ../../treadmill.bed --motif CGG -o ${i}_treadmill_out_DBSCAN/${i}.bin --store --plot --support ${support} --similarity 95.0 --threads 4 --affinity ${affinity}
				TREADMILL TRAP -i ${i}_treadmill_out_DBSCAN/${i}.bin -o ${i}_treadmill_out_DBSCAN

				predicted=$(bcftools query -f '%RALN,%AL1N' ${i}_treadmill_out_DBSCAN/TREADMILL.vcf.gz)
				rn=$(bcftools query -f '%RN\n' ${i}_treadmill_out_DBSCAN/TREADMILL.vcf.gz)
				en=$((${rn} + 30))
				expected="$rn,$en"
				echo -e ${i}"\t"${expected}"\t"${predicted}"\t"${l}X"\t"${j}%"\tDBSCAN" >> results_DBSCAN.txt

				TREADMILL RACE -fa ${genome} -bam ${i}_laser_out/sim.srt.bam -bed ../../treadmill.bed --motif CGG -o ${i}_treadmill_out_DENDOGRAM/${i}.bin --similarity 95.0 --threads 4 --hierarchical_clustering --dendogram --store
				python /home/davide/tools/TREADMILL/scripts/clusteranalysis.py -d ${i}_treadmill_out_DENDOGRAM/${i}.bin -s ${i}_treadmill_out_DENDOGRAM/simmatrix.bin -o ${i}_treadmill_out_SILHOUETTE
				n_clusters=$(tail -n+2 ${i}_treadmill_out_SILHOUETTE/silhouettescores.tsv | cut -f 1 | head -1)
				TREADMILL RACE -fa ${genome} -bam ${i}_laser_out/sim.srt.bam -bed ../../treadmill.bed --motif CGG -o ${i}_treadmill_out_AGGLOMERATIVE/${i}.bin --similarity 95.0 --threads 4 --hierarchical_clustering --clusters ${n_clusters} --store
				TREADMILL TRAP -i ${i}_treadmill_out_AGGLOMERATIVE/${i}.bin -o ${i}_treadmill_out_AGGLOMERATIVE

				predicted=$(bcftools query -f '%RALN,%AL1N' ${i}_treadmill_out_AGGLOMERATIVE/TREADMILL.vcf.gz)
				rn=$(bcftools query -f '%RN\n' ${i}_treadmill_out_AGGLOMERATIVE/TREADMILL.vcf.gz)
				en=$((${rn} + 30))
				expected="$rn,$en"
				echo -e ${i}"\t"${expected}"\t"${predicted}"\t"${l}X"\t"${j}%"\tAGGLOMERATIVE" >> results_AGGLOMERATIVE.txt

			done

		else

			echo "Skipping simulations"

		fi
		
		#calculate TPR for the 2 methods

		echo -e "TPR\terrors\tcoverage\tidentity\talgorithm\tstatus" > TPR_DBSCAN.txt
		echo -e "TPR\terrors\tcoverage\tidentity\talgorithm\tstatus" > TPR_AGGLOMERATIVE.txt

		er="0 1 2 3 4 5 6 7 8 9 10" #allow up to 10 errors

		for z in ${er}; do

			TP=$(tail -n+2 results_DBSCAN.txt | awk -v var=${z} '{split($3,a,","); if (a[2] != "." && a[2]>=50-var && a[2]<=50+var) print "0"; else print "1"}' | grep -c "0") #0 are correct; 1 incorrect
			total=$(tail -n+2 results_DBSCAN.txt | grep -c "^")
			TPR=$(echo ${TP}/${total} | bc -l)
			echo -e ${TPR}"\t"${z}"\t"${l}X"\t"${j}%"\tDBSCAN\tgreyzone" >> TPR_DBSCAN.txt

			TP=$(tail -n+2 results_AGGLOMERATIVE.txt | awk -v var=${z} '{split($3,a,","); if (a[2] != "." && a[2]>=50-var && a[2]<=50+var)  print "0"; else print "1"}' | grep -c "0") #0 are correct; 1 incorrect
			total=$(tail -n+2 results_AGGLOMERATIVE.txt | grep -c "^")
			TPR=$(echo ${TP}/${total} | bc -l)
			echo -e ${TPR}"\t"${z}"\t"${l}X"\t"${j}%"\tAGGLOMERATIVE\tgreyzone" >> TPR_AGGLOMERATIVE.txt
		
		done
		cd ..
	done 
	cd ..
done
cd ..

find .  -name "results*" -exec awk 'FNR>1 || NR==1' {} + > overview.tsv
find .  -name "TPR*" -exec awk 'FNR>1 || NR==1' {} + > tpr.tsv
