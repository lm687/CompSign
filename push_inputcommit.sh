#!/bin/bash

git add .
echo 'Everything added'
echo 'Add text for commit:'
read varname
git commit -m "$varname"
git push origin master
echo Done

