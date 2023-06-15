cat URL.txt | awk '{print "make -f ref-gen.mak -j4 PREFIX="$1,"URL="$2,">",$1".log 2>&1"}'
