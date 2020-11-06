#for i in Bnafr Brapa3.0 BolHDEM BniB Tpa Ath11
#do
#sed '1,5d' Bnafr.${i}.cds.psl|awk '{if($1/$11>=0.2){print "'$i'_"$10,$14}}'
#done|sed 's/\.[^ ]*//g'|sort|uniq|awk '{a[$1]++}END{for(i in a){print i,a[i]}}'|sed 's/_/ /'|awk '{if($3>1){s[$2]=1}else{s[$2]=s[$2]+0}if($1=="Bnafr"){a1[$2]=$3}if($1=="Brapa3.0"){a2[$2]=$3}if($1=="BolHDEM"){a3[$2]=$3}if($1=="BniB"){a4[$2]=$3}if($1=="Tpa"){a5[$2]=$3}if($1=="Ath11"){a6[$2]=$3}if($1=="Oglum"){a7[$2]=$3}if($1=="Omeri"){a8[$2]=$3}if($1=="Opunc"){a9[$2]=$3}if($1=="Obrac"){a10[$2]=$3}if($1=="Lperr"){a11[$2]=$3}}END{for(i in s){if(s[i]==0){print i,a1[i]+0,a2[i]+0,a3[i]+0,a4[i]+0,a5[i]+0,a6[i]+0,a7[i]+0,a8[i]+0,a9[i]+0,a10[i]+0,a11[i]+0}}}' >Bnafr.denovo.cds.status

for j in Bnafr Brapa3 BolHDEM BniB Tpa Ath11
do 
#sed '1,5d' Bnafr.${i}.dna.psl|sort -nk 1,1 |awk '{if($1/$11>=0.2){a[$10]++;b[$10]=$1/$11;c[$1]="'$i'_"$14}}END{for(i in a){print i,a[i],b[i],c[$1]}}'
sed '1,5d' Bnafr.${j}.dna.psl|sort -nk 1,1 |awk '{if($1/$11>=0.2){a[$10]++;b[$10]=$1/$11;c[$1]="'$j'"$14}}END{for(i in a){print i,a[i],b[i],c[$1]}}'
done |awk '{print $4"_"$1,$2,$3}' |awk '{if(a[$1]<$2){a[$1]=$2}if(b[$1]<$3){b[$1]=$3}}END{for(i in a){print i,a[i],b[i]}}' |sed 's/_/ /' |awk '{if($3>1){s[$2]=1}else{s[$2]=s[$2]+0}if($1=="BnafrchrC03"){a1[$2]=$4}if($1=="Brapa3A02"){a2[$2]=$4}if($1=="BolHDEMC5"){a3[$2]=$4}if($1=="BniBB5"){a4[$2]=$4}if($1=="Tpach3-1"){a5[$2]=$4}if($1=="Ath11Chr3"){a6[$2]=$4}}END{for(i in s){if(s[i]==0){print i,a1[i]+0,a2[i]+0,a3[i]+0,a4[i]+0,a5[i]+0,a6[i]+0}}}' >Bnafr.denovo.dna.status

#done #|sed 's/\./ /' |awk '{print $5"_"$1,$3,$4}' |awk '{if(a[$1]<$2){a[$1]=$2}if(b[$1]<$3){b[$1]=$3}}END{for(i in a){print i,a[i],b[i]}}' #|sed 's/_/ /'|awk '{if($3>1){s[$2]=1}else{s[$2]=s[$2]+0}if($1=="Bnafr"){a1[$2]=$4}if($1=="Brapa3"){a2[$2]=$4}if($1=="BolHDEM"){a3[$2]=$4}if($1=="BniB"){a4[$2]=$4}if($1=="Tpa"){a5[$2]=$4}if($1=="Ath11"){a6[$2]=$4}if($1=="Oglum"){a7[$2]=$4}if($1=="Omeri"){a8[$2]=$4}if($1=="Opunc"){a9[$2]=$4}if($1=="Obrac"){a10[$2]=$4}if($1=="Lperr"){a11[$2]=$4}}END{for(i in s){if(s[i]==0){print i,a1[i]+0,a2[i]+0,a3[i]+0,a4[i]+0,a5[i]+0,a6[i]+0,a7[i]+0,a8[i]+0,a9[i]+0,a10[i]+0,a11[i]+0}}}' >Bnafr.denovo.dna.status
