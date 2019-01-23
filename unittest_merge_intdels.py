#!/usr/bin/env python3

import unittest
import sys
sys.path.insert(0,'..')

from TEdist import merge


class Test_merge(unittest.TestCase):
    
    def setUp(self):
        pass
        
    def test_merge1(self):
        FT = {(100, 200) : [0.1, 0.1, 0.2], (102, 202) : [0.1, 0.3, 0.3]}
        ID = ((100, 200), [0.1, 0.1, 0.2])
        
        self.assertEqual(merge(freqtable=FT, intdel=ID), ( (101, 201), [0.2, 0.4, 0.5], [(100, 200), (102, 202)] ))
        
        
    def test_merge2(self):
        FT = {(100, 200) : [0.1, 0.1, 0.2], (104, 204) : [0.1, 0.3, 0.3]}
        ID = ((100, 200), [0.1, 0.1, 0.2])
        
        self.assertEqual(merge(freqtable=FT, intdel=ID), ( (100, 200), [0.1, 0.1, 0.2], [(100, 200)] ))
    
    
    def test_merge3(self):
        FT = {(100, 200) : [0.1, 0.1, 0.2], (102, 202) : [0.1, 0.3, 0.3], (98, 198) : [0.4, 0.4, 0.1]}
        ID = ((100, 200), [0.1, 0.1, 0.2])
        
        self.assertEqual(merge(freqtable=FT, intdel=ID), ( (100, 200), [0.6, 0.8, 0.6], [(100, 200), (102, 202), (98, 198)] ))
    
    
    def test_merge4(self):
        FT = {
            (100, 200) : [0.1, 0.1, 0.2],
            (102, 202) : [0.1, 0.3, 0.3],
            (98, 198) : [0.4, 0.4, 0.1],
            (200, 300) : [0.4, 0.4, 0.1]
            }
        ID = ((100, 200), [0.1, 0.1, 0.2])
        
        self.assertEqual(merge(freqtable=FT, intdel=ID), ( (100, 200), [0.6, 0.8, 0.6], [(100, 200), (102, 202), (98, 198)] ))



if __name__ == '__main__':
    unittest.main()
    
