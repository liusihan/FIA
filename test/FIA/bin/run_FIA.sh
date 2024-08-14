
nohup sh /home/hjb/FIA_github/test/FIA/bin/generate_evidence.0101.sh > /home/hjb/FIA_github/test/FIA/bin/generate_evidence.0101.log 2>&1 &
nohup sh /home/hjb/FIA_github/test/FIA/bin/generate_evidence.0102.sh > /home/hjb/FIA_github/test/FIA/bin/generate_evidence.0102.log 2>&1 &
nohup sh /home/hjb/FIA_github/test/FIA/bin/generate_evidence.0103.sh > /home/hjb/FIA_github/test/FIA/bin/generate_evidence.0103.log 2>&1 &
nohup sh /home/hjb/FIA_github/test/FIA/bin/generate_evidence.0104.sh > /home/hjb/FIA_github/test/FIA/bin/generate_evidence.0104.log 2>&1 &

for i in {1..900}
do
	record=`awk 'END{print NR}' /home/hjb/FIA_github/test/FIA/bin/record.txt`
	if [ "$record" -eq 5 ]
	then
		break
	else
		sleep 2m
	fi
done

nohup sh /home/hjb/FIA_github/test/FIA/bin/combine_nonsample_evidence.02.sh > /home/hjb/FIA_github/test/FIA/bin/combine_nonsample_evidence.02.log 2>&1 &

for i in {1..1000}
do
    record=`awk 'END{print NR}' /home/hjb/FIA_github/test/FIA/bin/record.txt`
    if [ "$record" -eq 6 ]
    then
        break
    else
        sleep 30s
    fi
done

nohup sh /home/hjb/FIA_github/test/FIA/bin/intervar_samplelist.sh > /home/hjb/FIA_github/test/FIA/bin/intervar_samplelist.log 2>&1 &


