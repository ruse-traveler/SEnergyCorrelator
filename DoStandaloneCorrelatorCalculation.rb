#!/usr/bin/env ruby
# -----------------------------------------------------------------------------
# 'DoStandaloneCorrelatorCalculation.rb'
# Derek Anderson
# 01.18.2024
#
# Short script to run the 'DoStandaloneCorrelatorCalculation.C' macro.
# -----------------------------------------------------------------------------

exec("root -b -q DoStandaloneCorrelatorCalculation.cxx")

# end -------------------------------------------------------------------------
