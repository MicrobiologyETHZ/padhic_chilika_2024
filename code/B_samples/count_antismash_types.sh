awk '/<div class=\"regbutton /{sub("\\s*<div class=\"regbutton ","",$0); sub(" r.*\">","",$0); split($0,line,","); for(i in line){counts[line[i]]+=1}}END{for(i in counts){print i,"\t",counts[i]}}' $1 | sort
