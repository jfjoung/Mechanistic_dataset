Acid_base = [
    # {Acidic form,Basic form, pKa, role}
    # Alykl Lithium
    {'Acid': '[#6;H1;+0:100].[Li;+1:99]','Base': '[#6;H0:100]-[Li:99]', 'pKa': 53, 'role' : ['B']},
    {'Acid': '[#6;H2;+0:100].[Li;+1:99]', 'Base': '[#6;H1:100]-[Li:99]', 'pKa': 51, 'role' : ['B']},
    {'Acid': '[#6;H3;+0:100].[Li;+1:99]', 'Base': '[#6;H2:100]-[Li:99]', 'pKa': 50, 'role' : ['B']},
    {'Acid': '[#6;H4;+0:100].[Li;+1:99]','Base': '[#6;H3:100]-[Li:99]', 'pKa': 48, 'role' : ['B']},
    
    {'Acid': '[#6;H1;+0:100]','Base': '[C;H0;D3;-1:100]', 'pKa': 53, 'role' : ['B']},
    {'Acid': '[#6;H2;+0:100]', 'Base': '[C;H1;D3;-1:100]', 'pKa': 51, 'role' : ['B']},
    {'Acid': '[#6;H3;+0:100]', 'Base': '[C;H2;D3;-1:100]', 'pKa': 50, 'role' : ['B']},
    {'Acid': '[#6;H4;+0:100]','Base': '[C;H3;D3;-1:100]', 'pKa': 48, 'role' : ['B']},
    
    # Lithium amide
    {'Acid': '[#7;H1;+0:100].[Li;+1:99]','Base': '[#7;H0:100]-[Li:99]', 'pKa': 37, 'role' : ['B']},
    {'Acid': '[#7;H2;+0:100].[Li;+1:99]','Base': '[#7;H1:100]-[Li:99]', 'pKa': 36, 'role' : ['B']},
    {'Acid': '[#7;H3;+0:100].[Li;+1:99]','Base': '[#7;H2:100]-[Li:99]', 'pKa': 35, 'role' : ['B']},
    {'Acid': '[#7;H1;+0:100]','Base': '[#7;H0;D2;-1:100]', 'pKa': 37, 'role' : ['B']},
    {'Acid': '[#7;H2;+0:100]','Base': '[#7;H1;D2;-1:100]', 'pKa': 36, 'role' : ['B']},
    {'Acid': '[#7;H3;+0:100]','Base': '[#7;H2;D2;-1:100]', 'pKa': 35, 'role' : ['B']},
    
    # Water, alcohol
    {'Acid': '[O;H3;+1:100]','Base':'[O;H2;+0:100]', 'pKa': -1.7, 'role' : ['A', 'B']},
    {'Acid': '[O;H2;+0:100]','Base':'[O;H1;-1:100]', 'pKa': 14, 'role' : ['A', 'B']},
    {'Acid': '[O;H2;+0:100].[Li,Na,K;+1:99]','Base':'[O;H1;+0:100]-[Li,Na,K;+0:99]', 'pKa': 14, 'role' : ['B']},
    {'Acid': '[C;!$([C]=[O]):100][O;H1;+0:99]', 'Base': '[C;!$([C]=[O]):100][O;H0;-1:99]', 'pKa':14, 'role' : ['A', 'B']},
    {'Acid': '[C;!$([C]=[O]):100][O;H1;+0:99].[Li,Na,K;+1:98]', 'Base': '[C;!$([C]=[O]):100][O;H0;+0:99]-[Li,Na,K;+0:98]', 'pKa':14, 'role' : ['B']},
    
    # Phenol
    {'Acid': '[c:100][O;H1;+0:99]', 'Base': '[c:100][O;H0;-1:99]', 'pKa':10, 'role' : ['A', 'B']},
    {'Acid': '[c:100][O;H1;+0:99].[Li,Na,K;+1:98]', 'Base': '[c:100][O;H0;+0:99]-[Li,Na,K;+0:98]', 'pKa':10, 'role' : ['B']},

    # thiol
    {'Acid': '[*:100][S;H1:99]', 'Base': '[*:100][S;H0;-1:99]', 'pKa':10, 'role' : ['A', 'B']},
    {'Acid': '[*:100][S;H1:99].[Li,Na,K;+1:98]', 'Base': '[*:100][S;H0;+0:99]-[Li,Na,K;+0:98]', 'pKa':10, 'role' : ['B']},
    
    #Selenol
    {'Acid': '[*:100][Se;H1:99]', 'Base': '[*:100][Se;H0;-1:99]', 'pKa':8, 'role' : ['A', 'B']},
    {'Acid': '[*:100][Se;H1:99].[Li,Na,K;+1:98]', 'Base': '[*:100][Se;H0;+0:99]-[Li,Na,K;+0:98]', 'pKa':8, 'role' : ['B']},
    
    # Sulfide
    {'Acid': '[S;H1;-1:100]','Base':'[S;H0;-2:100]', 'pKa': 12.89, 'role' : ['A', 'B']},
    {'Acid': '[Li,Na,K:100][S;H1;+0:99].[Li,Na,K;+1:98]','Base':'[Li,Na,K:100][S:99][Li,Na,K:98]', 'pKa': 12.89, 'role' : ['B']},
    
    {'Acid': '[S;H2;+0:100]','Base': '[S;H1;-1:100]', 'pKa': 7, 'role' : ['A', 'B']},
    {'Acid': '[Li,Na,K;+1:100].[S;H2;+0:99]','Base':'[Li,Na,K:100][S;H1;+0:99]', 'pKa': 7, 'role' : ['B']},
    
    #Phosphate
    {'Acid': '[O:100]=[P:99]([O;H0;-1:98])([O;H0;-1:97])[O;H1;+0:96]', 'Base':'[O:100]=[P:99]([O;H0;-1:98])([O;H0;-1:97])[O;H0;-1:96]', 'pKa': 12.32, 'role' : ['A', 'B']},
    {'Acid': '[O:100]=[P:99]([O;H0;-1:98])([O;H1;+0:97])[O;H1;+0:96]', 'Base': '[O:100]=[P:99]([O;H0;-1:98])([O;H0;-1:97])[O;H1;+0:96]', 'pKa': 7.21, 'role' : ['A', 'B']},
    {'Acid': '[O:100]=[P:99]([O;H1;+0:98])([O;H1;+0:97])[O;H1;+0:96]', 'Base': '[O:100]=[P:99]([O;H0;-1:98])([O;H1;+0:97])[O;H1;+0:96]', 'pKa': 2.12, 'role' : ['A', 'B']},
    
    {'Acid': '[O:100]=[P:99]([O;H1;+0:98])([O;H0;+0:97][*:95])[O;H0;+0:96][*:94]', 'Base': '[O:100]=[P:99]([O;H0;-1:98])([O;H0;+0:97][*:95])[O;H0;+0:96][*:94]', 'pKa': 1.5, 'role' : ['A', 'B']},
    {'Acid': '[O:100]=[P:99]([O;H1;+0:98])([O;H0;+0:97][*:95])[O;H1;+0:96]', 'Base': '[O:100]=[P:99]([O;H0;-1:98])([O;H0;+0:97][*:95])[O;H1;+0:96]', 'pKa': 1.0, 'role' : ['A', 'B']},
        
    #Cyanate
    {'Acid': '[C;H1;+0:100]#[N;H0;+0:99]', 'Base': '[C;H0;-1:100]#[N;H0;+0:99]', 'pKa': 9.40, 'role' : ['A', 'B']},
    {'Acid': '[C;H1;+0:100]#[N;H0;+0:99].[Li,Na,K;+1:98]', 'Base': '[Li,Na,K:98][C;H0;+0:100]#[N;H0;+0:99]', 'pKa': 9.40, 'role' : ['B']},
    
    # Ammonia, amine
    {'Acid': '[N;D0;H4;+1:100]' , 'Base': '[N;D0;H3;+0:100]', 'pKa': 9.25, 'role' : ['A', 'B']},
    {'Acid' :'[N;D0;H4;+0:100][Cl,Br,I:99]' , 'Base': '[N;D0;H3;+0:100].[Cl,Br,I;-1:99]', 'pKa': 9.25, 'role' : ['A']},
    {'Acid' :'[C;!$(C=*):100]-[N;D3;H1;+1:99](-[C;!$(C=*):98])-[C;!$(C=*):97]' , 'Base': '[C;!$(C=*):100]-[N;D3;H0;+0:99](-[C;!$(C=*):98])-[C;!$(C=*):97]', 'pKa': 10.8, 'role' : ['A', 'B']},
    {'Acid' :'[C;!$(C=*):100]-[N;D2;H2;+1:99]-[C;!$(C=*):98]', 'Base': '[C;!$(C=*):100]-[N;D2;H1;+0:99]-[C;!$(C=*):98]', 'pKa': 11, 'role' : ['A', 'B']},
    {'Acid' :'[C;!$(C=*):100]-[N;D1;H3;+1:99]' , 'Base': '[C;!$(C=*):100]-[N;D1;H2;+0:99]', 'pKa': 10.6, 'role' : ['A', 'B']},
    
    {'Acid' :'[O:100]=[C;R:99]-[N;H1;+0;R:98]-[C;R:97]=[O:96]' , 'Base': '[O:100]=[C;R:99]-[N;H0;-1;R:98]-[C;R:97]=[O:96]', 'pKa': 9.6, 'role' : ['A', 'B']},    

    # Amide
    {'Acid' :'[N;H1;+1:100][C:98]=[O:97]' , 'Base': '[N;H0;+0:100][C:98]=[O:97]', 'pKa': -0.3, 'role' : ['A', 'B']},
    
    
    #Schiff base
    {'Acid': '[C;+0:100]=[N;H1;+1:99][C;+0:98]', 'Base': '[C;+0:100]=[N;+0;H0:99][C;+0:98]', 'pKa': 14, 'role' : ['A', 'B']},
    
    # Heterocyclic amine
    #Pyrimidine
    {'Acid': '[n&!$([#7]~[#7]);+1;H1:100]1[c:99][c:98][c:97][n:96][c:95]1', 'Base': '[n&!$([#7]~[#7]);H0:100]1[c:99][c:98][c:97][n:96][c:95]1', 'pKa': 1.3, 'role' : ['A', 'B']},
    #Pyridine
    {'Acid': '[n&!$([#7]~[#7]);+1;H1:100]1[c:99][c:98][c:97][c:96][c:95]1', 'Base': '[n&!$([#7]~[#7]);H0:100]1[c:99][c:98][c:97][c:96][c:95]1', 'pKa': 5.25, 'role' : ['A', 'B']},
    #Imidazole
    {'Acid': '[c:100]1[n;H1:99][c:96][c:98][n&!$([#7]~[#7]);H1;+1:97]1', 'Base': '[c:100]1[n;H1:99][c:96][c:98][n&!$([#7]~[#7]);H0;+0:97]1', 'pKa': 6.29, 'role' : ['A', 'B']},
    #Pyrrole
    {'Acid': '[c:100]1[n&!$([#7]~[#7]);H2;+1:99][c:96][c:98][c:97]1', 'Base': '[c:100]1[n&!$([#7]~[#7]);H1;+0:99][c:96][c:98][c:97]1', 'pKa': 0.4, 'role' : ['A', 'B']},
    # DBU, 1,8~Diazabicyclo(5.4.0)undec~7~ene
    {'Acid': '[C:100][N&!$([#7]~[#7]);H1;+1:99]=[C;H0:98][N;H0:97]','Base':'[C:100][N&!$([#7]~[#7]);+0;H0:99]=[C;H0:98][N;H0:97]', 'pKa': 13.4, 'role' : ['A', 'B']},
    
    #Aniline
    {'Acid': '[N&!$([#7]~[#7]);H3;+1:100][c:99]1[c:98][c:97][c:96][c:95][c:94]1', 'Base': '[N&!$([#7]~[#7]);H2;+0:100][c:99]1[c:98][c:97][c:96][c:95][c:94]1', 'pKa': 4.63, 'role' : ['A', 'B']},
    {'Acid': '[*:101][N&!$([#7]~[#7]);H2;+1:100][c:99]1[c:98][c:97][c:96][c:95][c:94]1', 'Base': '[*:101][N&!$([#7]~[#7]);H1;+0:100][c:99]1[c:98][c:97][c:96][c:95][c:94]1', 'pKa': 4.63, 'role' : ['A', 'B']},
    {'Acid': '[*:101][N&!$([#7]~[#7]);H1;+1:100]([*:102])[c:99]1[c:98][c:97][c:96][c:95][c:94]1', 'Base': '[*:101][N&!$([#7]~[#7]);H0;+0:100]([*:102])[c:99]1[c:98][c:97][c:96][c:95][c:94]1', 'pKa': 4.63, 'role' : ['A', 'B']},
    
    #Inorganic acid and base
    # hydrocyanic acid
    {'Acid': '[C;H1;+0:100]#[N:99]' , 'Base':'[C;H0;-1:100]#[N:99]' , 'pKa': 9.40, 'role' : ['A', 'B']},
    {'Acid': '[C;H1;+0:100]#[N:99].[Li,Na,K;+1:98]' , 'Base':'[Li,Na,K;+0:98][C;H0;+0:100]#[N:99]' , 'pKa': 9.40, 'role' : ['B']},
    # Hypochlorous acid
    {'Acid': '[O;H1;+0:100][Cl;+0:99]', 'Base': '[O;H0;-1:100][Cl;+0:99]', 'pKa': 7.46, 'role' : ['A', 'B']},
    {'Acid': '[O;H1;+0:100][Cl;+0:99].[Li,Na,K;+1:98]', 'Base': '[Li,Na,K;+0:98][O;H0;+0:100][Cl;+0:99]', 'pKa': 7.46, 'role' : ['B']},
    #Sulfuric acid
    {'Acid': '[O:100]=[S:99]([O;H1;+0:98])([O;H1;+0:97])=[O:96]', 'Base': '[O:100]=[S:99]([O;H1;+0:98])([O;H0;-1:97])=[O:96]', 'pKa': -10, 'role' : ['A', 'B']},
    {'Acid': '[O:100]=[S:99]([O;H1;+0:98])([O;H1;+0:97])=[O:96].[Li,Na,K;+1:95]', 'Base': '[O:100]=[S:99]([O;H1;+0:98])([O;H0;+0:97][Li,Na,K;+0:95])=[O:96]', 'pKa': -10, 'role' : ['B']},
    {'Acid': '[O:100]=[S:99]([O;H1;+0:98])([O;H0;-1:97])=[O:96]', 'Base': '[O:100]=[S:99]([O;H0;-1:98])([O;H0;-1:97])=[O:96]', 'pKa': 1.92, 'role' : ['A', 'B']},
    {'Acid': '[O:100]=[S:99]([O;H0;+0:98][Li,Na,K;+0:94])([O;H1;+0:97])=[O:96].[Li,Na,K;+1:95]', 'Base': '[O:100]=[S:99]([O;H0;+0:98][Li,Na,K;+0:94])([O;H0;+0:97][Li,Na,K;+0:95])=[O:96]', 'pKa': 1.92, 'role' : ['B']},
    {'Acid': '[O:100]=[S:99]([O;H0;-1:98])([O;H1;+0:97])=[O:96].[Li,Na,K;+1:95]', 'Base': '[O:100]=[S:99]([O;H0;-1:98])([O;H0;+0:97][Li,Na,K;+0:95])=[O:96]', 'pKa': 1.92, 'role' : ['B']},
    
    #sulfurous acid
    {'Acid': '[O:100]=[S;D3:99]([O;H1;+0:98])([O;H1;+0:97])', 'Base': '[O:100]=[S;D3:99]([O;H1;+0:98])([O;H0;-1:97])', 'pKa': 1.82, 'role' : ['A', 'B']},
    {'Acid': '[O:100]=[S;D3:99]([O;H1;+0:98])([O;H1;+0:97]).[Li,Na,K;+1:95]', 'Base': '[O:100]=[S;D3:99]([O;H1;+0:98])([O;H0;+0:97][Li,Na,K;+0:95])', 'pKa': 1.82, 'role' : ['B']},
    {'Acid': '[O:100]=[S;D3:99]([O;H1;+0:98])([O;H0;-1:97])', 'Base': '[O:100]=[S;D3:99]([O;H0;-1:98])([O;H0;-1:97])', 'pKa': 7, 'role' : ['A', 'B']},
    {'Acid': '[O:100]=[S;D3:99]([O;H0;+0:98][Li,Na,K;+0:94])([O;H1;+0:97]).[Li,Na,K;+1:95]', 'Base': '[O:100]=[S;D3:99]([O;H0;+0:98][Li,Na,K;+0:94])([O;H0;+0:97][Li,Na,K;+0:95])', 'pKa': 7, 'role' : ['B']},
    {'Acid': '[O:100]=[S;D3:99]([O;H0;-1:98])([O;H1;+0:97]).[Li,Na,K;+1:95]', 'Base': '[O:100]=[S;D3:99]([O;H0;-1:98])([O;H0;+0:97][Li,Na,K;+0:95])', 'pKa': 7, 'role' : ['B']},
    
    #organic sulfonic acid, TsOH
    {'Acid':'[O:100]=[S:99](-[#6:98])(-[O;H1;+0:97])=[O:96]' , 'Base': '[O:100]=[S:99](-[#6:98])([O;H0;-1:97])=[O:96]', 'pKa': -2.8, 'role' : ['A', 'B']},
    {'Acid':'[O:100]=[S:99](-[#6:98])(-[O;H1;+0:97])=[O:96].[Li,Na,K;+1:95]' , 'Base': '[O:100]=[S:99](-[#6:98])([O;H0;+0:97][Li,Na,K;+0:95])=[O:96]', 'pKa': -2.8, 'role' : ['B']},
    
    #Carboxylate
    {'Acid': '[C:100](=[O:99])[O;H1;+0:98]' , 'Base': '[C:100](=[O:99])[O;H0;-1:98]', 'pKa': 4.74, 'role' : ['A', 'B']},
    {'Acid': '[C:100](=[O:99])[O;H1;+0:98].[Li,Na,K;+1:97]' , 'Base': '[C:100](=[O:99])[O;H0;+0:98][Li,Na,K;+0:97]', 'pKa': 4.74, 'role' : ['B']},
    {'Acid': '[O;H1;+0:100][C:99]([C:98]([F:97])([F:96])[F:95])=[O:94]' , 'Base': '[O;H0;-1:100][C:99]([C:98]([F:97])([F:96])[F:95])=[O:94]', 'pKa': -0.3, 'role' : ['A', 'B']},
    {'Acid': '[O;H1;+0:100][C:99]([C:98]([F:97])([F:96])[F:95])=[O:94].[Li,Na,K;+1:93]' , 'Base': '[Li,Na,K;0:93][O;H0;+0:100][C:99]([C:98]([F:97])([F:96])[F:95])=[O:94]', 'pKa': -0.3, 'role' : ['B']},

    
    #hydrazoic acid
    {'Acid': '[N;H1;+0:100]=[N;+1:99]=[N;-1:98]', 'Base': '[N;H0;-1:100]=[N;+1:99]=[N;-1:98]', 'pKa': 4.72, 'role' : ['A', 'B']},
    {'Acid': '[N;H1;+0:100]=[N;+1:99]=[N;-1:98].[Li,Na,K;+1:97]', 'Base': '[Li,Na,K;+0:97][N;H0;+0:100]=[N;+1:99]=[N;-1:98]', 'pKa':4.72, 'role' : ['B']},
    
    #Nitric acid
    {'Acid': '[O:100]=[N;+1:99]([O;-1:98])[O;H1;+0:97]', 'Base': '[O:100]=[N;+1:99]([O;-1:98])[O;H0;-1:97]', 'pKa': -1.4, 'role' : ['A', 'B']},
    
    #Nitrous acid
    {'Acid': '[O:100]=[N;+1:99][O;H1;+0:98]', 'Base': '[O:100]=[N;+1:99][O;H0;-1:98]', 'pKa': 3.4, 'role' : ['A', 'B']},
    {'Acid': '[O:100]=[N;+1:99][O;H1;+0:98].[Li,Na,K;+1:97]', 'Base': '[O:100]=[N;+1:99][O;H0;+0:98][Li,Na,K;+0:97]', 'pKa': 3.4, 'role' : ['B']},
    
    #Perchloric acid
    {'Acid': '[O;H1;+0:100][Cl:99](=[O:98])(=[O:97])=[O:96]', 'Base': '[O;H0;-1:100][Cl:99](=[O:98])(=[O:97])=[O:96]', 'pKa': -15.2 , 'role' : ['A', 'B']},
    
    #hydrohalic acid
    {'Acid': '[Cl,Br,I;H1;+0:100]' , 'Base': '[Cl,Br,I;H0;-1:100]', 'pKa': -6, 'role' : ['A', 'B']},
    {'Acid': '[F;H1;+0:100]' , 'Base': '[F;H0;-1:100]', 'pKa': 3.14, 'role' : ['A', 'B']},
    {'Acid': '[F;H1;+0:100].[Li,Na,K,Cs;+1:95]' , 'Base': '[F;H0;+0:100][Li,Na,K,Cs;+0:95]', 'pKa': 3.14, 'role' : ['B']},
    
    #Iodic acid
    {'Acid': '[O;H1;+0:100][I:99](=[O:98])=[O:97]', 'Base': '[O;H0;-1:100][I:99](=[O:98])=[O:97]', 'pKa': 0.77, 'role' : ['A', 'B']},
    
    
    #Carbonate
    {'Acid':'[O:100]=[C:99]([O;H1;+0:98])[O;H0;-1:97]' , 'Base': '[O:100]=[C:99]([O;H0;-1:98])[O;H0;-1:97]' , 'pKa': 10.33, 'role' : ['A', 'B']},
    {'Acid':'[O:100]=[C:99]([O;H1;+0:98])[O;H0;-1:97].[Li,Na,K,Cs;+1:95]' , 'Base': '[O:100]=[C:99]([O;H0;+0:98][Li,Na,K,Cs;+0:95])[O;H0;-1:97]' , 'pKa': 10.33, 'role' : ['B']},
    {'Acid':'[O:100]=[C:99]([O;H1;+0:98])[O;H0;+0:97][Li,Na,K,Cs;+0:94].[Li,Na,K,Cs;+1:95]' , 'Base': '[O:100]=[C:99]([O;H0;+0:98][Li,Na,K,Cs;+0:95])[O;H0;+0:97][Li,Na,K,Cs;+0:94]' , 'pKa': 10.33, 'role' : ['B']},
    {'Acid':'[O:100]=[C:99]([O;H1;+0:98])[O;H1;+0:97]' , 'Base': '[O:100]=[C:99]([O;H1;+0:98])[O;H0;-1:97]' , 'pKa': 6.37, 'role' : ['A', 'B']},
    {'Acid':'[O:100]=[C:99]([O;H1;+0:98])[O;H1;+0:97].[Li,Na,K,Cs;+1:94]' , 'Base': '[O:100]=[C:99]([O;H1;+0:98])[O;H0;+0:97][Li,Na,K,Cs;+0:94]' , 'pKa': 6.37, 'role' : ['B']},
 
    
    
    ]

   
    
#     {'Acid': , 'Base': , 'pKa': , 'role' : ['A', 'B']},
#     {'Acid': , 'Base': , 'pKa': , 'role' : ['A', 'B']},
#     {'Acid': , 'Base': , 'pKa': , 'role' : ['A', 'B']},
#     {'Acid': , 'Base': , 'pKa': , 'role' : ['A', 'B']},
#     {'Acid': , 'Base': , 'pKa': , 'role' : ['A', 'B']},
#     {'Acid': , 'Base': , 'pKa': , 'role' : ['A', 'B']},
#     {'Acid': , 'Base': , 'pKa': , 'role' : ['A', 'B']},
#     {'Acid': , 'Base': , 'pKa': , 'role' : ['A', 'B']},