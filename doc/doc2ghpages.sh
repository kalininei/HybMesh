#!/usr/bin/bash

# publishes doc/build/html to gh-pages.
# should be called from project root directory.

# build docs
cd doc
make
cd ..

# copy html directory outside main repository to a tmp destination.
cp -R doc/build/html ../_ghpages_temporary_directory

# switch to gh-pages branch
git stash
git checkout gh-pages
rm -rf *
