#!/bin/bash

# $1 : fastafile

outaligned=${1%.*}"_aligned.fa"
if [ ! -f $outaligned ]
    then /home/adrien/source/muscle3.8.425_i86linux64 -in $1 -out $outaligned || exit 1
fi

outcleaned=${outaligned%.*}"_cleaned.fa"
if [ ! -f $outcleaned ]
    then Aln-remove-gaps -m 30 $outaligned && mv $outaligned".cleaned" $outcleaned || exit 1
fi

outrenamed=${outaligned%.*}"_renamed.fa"
if [ ! -f $outrenamed ]
    then awk '/^>/{print ">a" ++i; next}{print}' < $outcleaned > $outrenamed || exit 1
fi

outrenamed_phyl=${outrenamed%.*}".phylip"
if [ ! -f $outrenamed_phyl ]
    then fasta2phylip.py $outrenamed $outrenamed_phyl || exit 1
fi

outcleaned_phyl=${outcleaned%.*}".phylip"
if [ ! -f $outcleaned_phyl ]
    then fasta2phylip.py $outcleaned $outcleaned_phyl || exit 1
fi

fastree=${outrenamed%.*}".nwk"
if [ ! -f $fastree ]
    then FastTree $outrenamed_phyl > $fastree || exit 1
fi

model=${1%.*}".model"
if [ ! -f $model ]
    then java -jar /home/adrien/source/prottest-3.4-20140123/prottest-3.4.jar -i $outrenamed_phyl -t $fastree -S 0 -all-distributions -F -AIC -BIC > $model 
fi

m=$(grep -oP "(?<=Best model according to AIC: ).*" $model)
arr=(${m//+/ })

mod=${arr[0]}
unset arr[0]
plus=" "
for i in "${arr[@]}"; do if [ "$i" == "I" ]; then plus+="-v e "; elif [ "$i" == "G" ]; then plus+="-a e "; elif [ "$i" == "F" ]; then plus+="-f e "; fi done

phytree=$outcleaned_phyl"_phyml_tree.txt"
phyout=${phytree%.*}".nwk"
if [ ! -f $phytree ]
    then phyml -i "$outcleaned_phyl" -d "aa" -m "$mod" "$plus" -b -1 && mv $phytree $phyout || exit 1
fi


