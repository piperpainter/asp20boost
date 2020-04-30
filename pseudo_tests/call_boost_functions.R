### --------------------------------------------------------------
# Function-Calls -------------------------------------------------

# This file simply tries to asses the boost-function we intend to push.


# Load our own package -------------------------------------------
library(asp20boost)


# Construct pseudo-model to pass to boost-function ---------------
x <- 1

# Call Functions -------------------------------------------------
boost_levin(x)
try(boost_sebastian(x))
try(boost_johannes(x))