Overview of the toy example dataset for testing:
================================================
                    10        20        30        40    
Coord:    01234567890123456789012345678901234567890
Domain    <---------->   <-------------------->        (0,12), (15,37)
Boundary              ###                              (12,15)
GeneC,B,D,A -->     --->    -->     ----->             C:(2,5), B:(10,14), D:(18,20), A:(26,32)
Enhacer     **                      **                 (2,4), (26,28) 
cnv1               ==========                          (9,19)
cnv2              =========================            (8,33)
cnv3                  =======                          (12,19)
cnv4                    ==                             (14,16)
                    10        20        30        40    
          01234567890123456789012345678901234567890


Example phenotype ontology (EPO) data for testing:

                    EP00
                   /    \       
                EP01     EP02        
               /       /    \       
            EP03    EP04    EP05
                       \    /  \
                        EP06    EP07

Term    p   IC=-log(p)
========================
EP01    1   0
EP02    .25 1.39
EP03    .75 0.29
EP04    .25 1.39
EP05    .75 0.29
EP06    0   1.39 (Inf)
EP07    .25 1.39

CNV phenotypes: EP06
CNV targetTerm: EP06

Gene    Annotation  PhenoMatchScore (to EP:06)
==================================
GeneA   EP04,EP05   1.68
GeneB   EP07        0.29
GeneC   EP03        0
GeneD   EP05        0.29

