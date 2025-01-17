#!/bin/sh 

for file_name in *.gjf; do
    awk '/[0-9]\.[0-9][0-9]/ { a[++n] = $0 }
      END { print n; print; 
        for(i=1; i<=n; ++i) print a[i] }' "$file_name" > "${file_name%.gjf}.xyz"
done