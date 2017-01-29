"""
This module contains functions to calculate a running mean, median and mode.
"""

#  The MIT License (MIT)
#
#  Copyright (c) 2013 Virgil Stokes
#  http://code.activestate.com/recipes/578480-running-median-mean-and-mode/
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
#
#  NOTE: The MIT License is limited to this document

######################### Python 2 and 3 compatibility #########################
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
from future.utils import viewitems, viewkeys, viewvalues, listitems, listvalues

try:
    #python 2.7
    from itertools import izip as zip
except ImportError:
    #python 3 series
    pass
################################################################################

from bisect import insort, bisect_left
from collections import deque, Counter
from itertools import islice


def runningMedian(seq, M):
    """
     Purpose: Find the median for the points in a sliding window (odd number in size)
              as it is moved from left to right by one point at a time.
      Inputs:
            seq -- list containing items for which a running median (in a sliding window)
                   is to be calculated
              M -- number of items in window (window size) -- must be an integer > 1
      Otputs:
         medians -- list of medians with size N - M + 1
       Note:
         1. The median of a finite list of numbers is the "center" value when this list
            is sorted in ascending order.
         2. If M is an even number the two elements in the window that
            are close to the center are averaged to give the median (this
            is not by definition)
    """
    seq = iter(seq)
    s = []
    m = M // 2 #// does a truncated division like integer division in Python 2

    # Set up list s (to be sorted) and load deque with first window of seq
    s = [item for item in islice(seq,M)]
    d = deque(s)

    # Simple lambda function to handle even/odd window sizes
    median = lambda : s[m] if bool(M&1) else (s[m-1]+s[m])*0.5

    # Sort it in increasing order and extract the median ("center" of the sorted window)
    s.sort()
    medians = [median()]

    # Now slide the window by one point to the right for each new position (each pass through
    # the loop). Stop when the item in the right end of the deque contains the last item in seq
    for item in seq:
        old = d.popleft()          # pop oldest from left
        d.append(item)             # push newest in from right
        del s[bisect_left(s, old)] # locate insertion point and then remove old
        insort(s, item)            # insert newest such that new sort is not required
        medians.append(median())
    return medians


def runningMean(seq, N, M):
    """
     Purpose: Find the mean for the points in a sliding window (fixed size)
              as it is moved from left to right by one point at a time.
      Inputs:
          seq -- list containing items for which a mean (in a sliding window) is
                 to be calculated (N items)
            N -- length of sequence
            M -- number of items in sliding window
      Otputs:
        means -- list of means with size N - M + 1
    """
    # Load deque (d) with first window of seq
    d = deque(seq[0:M])
    means = [sum(d) / len(d)] # contains mean of first window
    # Now slide the window by one point to the right for each new position (each pass through
    # the loop). Stop when the item in the right end of the deque contains the last item in seq
    for item in islice(seq, M, N):
        old = d.popleft()            # pop oldest from left
        d.append(item)               # push newest in from right
        m = sum(d) / len(d)
        means.append(m)     # mean for current window
    return means
