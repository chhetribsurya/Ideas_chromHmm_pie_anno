In [85]: df[df.Mnemonics.str.match("\w+Ctcf*")]
Out[85]: 
   Mnemonics                                         Rationale State_index
0    TssCtcf               CTCF at transcription start sites\n           0
27  Gen3Ctcf  CTCF enriched in transcription 3’ end of genes\n          27
33  PromCtcf           CTCF enriched in regions flanking TSS\n          33

In [86]: Ctcf_df = pd.concat([df[df.Mnemonics.str.match("\w+Ctcf*")],df[df.Mnemonics.str.match("\w?Ctcf*")]], axis = 0)

In [87]: Ctcf_df
Out[87]: 
   Mnemonics                                          Rationale State_index
0    TssCtcf                CTCF at transcription start sites\n           0
27  Gen3Ctcf   CTCF enriched in transcription 3’ end of genes\n          27
33  PromCtcf            CTCF enriched in regions flanking TSS\n          33
23      Ctcf  Distal CTCF/candidate insulator without open c...          23
34     CtcfO  Distal CTCF/candidate insulator with open chro...          34

In [88]: df[df.Mnemonics.str.match("\w+Repr*")]
Out[88]: 
   Mnemonics                              Rationale State_index
31  LowReprW  Low signal/weak polycomb repression\n          31

In [90]: df[df.Mnemonics.str.match("\w+Repr*|\w?Repr*")]
Out[90]: 
   Mnemonics                                    Rationale State_index
4      Repr1                 Strong polycomb repression\n           4
10     Repr2                 Strong polycomb repression\n          10
21     ReprD  Polycomb repression with Duke DNase sites\n          21
31  LowReprW        Low signal/weak polycomb repression\n          31

In [91]: Polycomb_df = df[df.Mnemonics.str.match("\w+Repr*|\w?Repr*")]

In [92]: df[df.Mnemonics.str.match("\w+Dnase*|\w+Faire")]
Out[92]: 
Empty DataFrame
Columns: [Mnemonics, Rationale, State_index]
Index: []

In [93]: df[df.Mnemonics.str.match("\w+Dnase*|\w+Faire*")]
Out[93]: 
Empty DataFrame
Columns: [Mnemonics, Rationale, State_index]
Index: []

In [94]: df[df.Mnemonics.str.match("\w+Dnase*")]
Out[94]: 
Empty DataFrame
Columns: [Mnemonics, Rationale, State_index]
Index: []

In [95]: df.Mnemonics
Out[95]: 
0      TssCtcf
1      FaireW2
2      FaireW1
3        H4K20
4        Repr1
5       PromF1
6         Elon
7       EnhWF3
8      DnaseD1
9        Quies
10       Repr2
11        Low2
12     DnaseD2
13        Gen5
14      EnhWF2
15      EnhWF1
16         Tss
17        Zero
18      PromF2
19        Gen3
20        Pol2
21       ReprD
22        TssF
23        Ctcf
24       PromP
25        Art2
26        TssW
27    Gen3Ctcf
28        Art1
29        EnhF
30        Low1
31    LowReprW
32         Enh
33    PromCtcf
34       CtcfO
35        EnhW
36       Table
Name: Mnemonics, dtype: object

In [96]: df[df.Mnemonics.str.match("Dnase*|Faire*")]
Out[96]: 
   Mnemonics                           Rationale State_index
1    FaireW2  Modest Faire/Control enrichments\n           1
2    FaireW1  Modest Faire/Control enrichments\n           2
8    DnaseD1              Primarily Duke DNase\n           8
12   DnaseD2              Primarily Duke DNase\n          12

In [97]: Assay_faireDnase_df = df[df.Mnemonics.str.match("Dnase*|Faire*")]

In [98]: df[df.Mnemonics.str.match("Art* | Quies*")]
Out[98]: 
Empty DataFrame
Columns: [Mnemonics, Rationale, State_index]
Index: []

In [99]: df[df.Mnemonics.str.match("Art*|Quies*")]
Out[99]: 
   Mnemonics                                Rationale State_index
9      Quies                        Heterochromatin\n           9
25      Art2  Potential CNV or repetitive artifacts\n          25
28      Art1  Potential CNV or repetitive artifacts\n          28

In [100]: Repeats_df = df[df.Mnemonics.str.match("Art*|Quies*")]

In [101]: Zero_dead_df = df[df.Mnemonics.str.match("Zero*")]

In [102]: Zero_dead_df
Out[102]: 
   Mnemonics    Rationale State_index
17      Zero  Dead zone\n          17


In [121]: df[df.Mnemonics.str.match("Tss[\w]?")]
Out[121]: 
   Mnemonics                                    Rationale State_index
0    TssCtcf          CTCF at transcription start sites\n           0
16       Tss  Transcription start sites/active promoter\n          16
22      TssF              Active promoter, flanking TSS\n          22
26      TssW                         Weak promoter, TSS\n          26

In [122]: df[df.Mnemonics.str.match("Tss[\w]?$")]
Out[122]: 
   Mnemonics                                    Rationale State_index
16       Tss  Transcription start sites/active promoter\n          16
22      TssF              Active promoter, flanking TSS\n          22
26      TssW                         Weak promoter, TSS\n          26

