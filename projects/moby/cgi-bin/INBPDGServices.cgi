#!/bin/sh

# First we need to load environmental variables in order to use perl
. /etc/profile

# Then, go!
exec perl "${0}.pl" "$@"
