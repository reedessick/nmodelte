python_alias="`which python | tr " " "\n" | tail -n 1`"
for filename in `ls bin/*py`
	do 
	sed -i "s;python_alias;${python_alias};" ${filename}
	done
