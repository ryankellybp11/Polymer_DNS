  Ð7  ]   k820309    Ö          24.0        ÐRf                                                                                                          
       code/src/dns/derivatives.f90 DERIVS #         @                                                               #DERIVS!GRADV!WAVES    #U    #V    #W    #U11 	   #U12 
   #U13    #U21    #U22    #U23    #U31    #U32    #U33    #LU    #LV    #LW                                                                                                                                                         #GRADV%WAVES%WAVX    #GRADV%WAVES%WAVZ    #GRADV%WAVES%C                                                                                
      p          p            p                                                                                                             
      p          p            p                                                                                                             
      p          p            p                                            D @                                              @                  p A        p          p          p            p          p          p                                    D @                                              @                  p A        p          p          p            p          p          p                                    D @                                              @                  p A        p          p          p            p          p          p                                    D                                          	      @                  p A        p          p          p            p          p          p                                    D @                                        
      @                  p A        p          p          p            p          p          p                                    D                                                @              	    p A        p          p          p            p          p          p                                    D                                                @              
    p A        p          p          p            p          p          p                                    D @                                              @                  p A        p          p          p            p          p          p                                    D                                                @                  p A        p          p          p            p          p          p                                    D                                                @                  p A        p          p          p            p          p          p                                    D @                                              @                  p A        p          p          p            p          p          p                                    D                                                @                  p A        p          p          p            p          p          p                                    D                                                @                  p A        p          p          p            p          p          p                                    D                                                @                  p A        p          p          p            p          p          p                                    D                                                @                  p A        p          p          p            p          p          p                          #         @                                                               #DERIVS!CDERIV!SOLVER    #F    #DF                                                                                                                                          0                     #CDERIV%SOLVER%GAIN    #CDERIV%SOLVER%UGAIN    #CDERIV%SOLVER%THETA    #CDERIV%SOLVER%ALPHA    #CDERIV%SOLVER%BETA    #CDERIV%SOLVER%DYDE                                                                  
                                                                    
                                                                    
                                                                    
                                                                     
                                                             (       
                                                                 @                  p A        p          p          p            p          p          p                                    D                                                @                  p A        p          p          p            p          p          p                          #         @                                                                #WRKC     #WRK1 !   #IY "                                                                                                                                                 @                  p A        p          p          p            p          p          p                                    D                                          !      @                  p A        p          p          p            p          p          p                                                                              "            #         @                                            #                   #DERIVS!C1DERBW!SOLVER $   #WRKC +   #WRK1 ,   #IY -                                                                                                                                         $     0                     #C1DERBW%SOLVER%GAIN %   #C1DERBW%SOLVER%UGAIN &   #C1DERBW%SOLVER%THETA '   #C1DERBW%SOLVER%ALPHA (   #C1DERBW%SOLVER%BETA )   #C1DERBW%SOLVER%DYDE *                                                    %             
                                                        &            
                                                        '            
                                                        (            
                                                        )             
                                                        *     (       
                                                           +      @                  p A        p          p          p            p          p          p                                    D                                          ,      @                  p A        p          p          p            p          p          p                                                                              -            #         @                                            .                   #DERIVS!C2DERBW!SOLVER /   #WRKC 6   #WRK1 7   #IY 8                                                                                                                                         /     0                     #C2DERBW%SOLVER%GAIN 0   #C2DERBW%SOLVER%UGAIN 1   #C2DERBW%SOLVER%THETA 2   #C2DERBW%SOLVER%ALPHA 3   #C2DERBW%SOLVER%BETA 4   #C2DERBW%SOLVER%DYDE 5                                                    0             
                                                        1            
                                                        2            
                                                        3            
                                                        4             
                                                        5     (       
                                                           6      @                  p A        p          p          p            p          p          p                                    D                                          7      @                  p A        p          p          p            p          p          p                                                                              8            #         @                                            9                    #WRKC :   #WRK1 ;   #IY <                                                                                                                                          :      @                  p A        p          p          p            p          p          p                                    D                                          ;      @                  p A        p          p          p            p          p          p                                                                              <            #         @                                            =                   #DERIVS!C1DERTW!SOLVER >   #WRKC E   #WRK1 F   #IY G                                                                                                                                         >     0                     #C1DERTW%SOLVER%GAIN ?   #C1DERTW%SOLVER%UGAIN @   #C1DERTW%SOLVER%THETA A   #C1DERTW%SOLVER%ALPHA B   #C1DERTW%SOLVER%BETA C   #C1DERTW%SOLVER%DYDE D                                                    ?             
                                                        @            
                                                        A            
                                                        B            
                                                        C             
                                                        D     (       
                                                           E      @                   p A        p          p          p            p          p          p                                    D                                          F      @              !    p A        p          p          p            p          p          p                                                                              G            #         @                                            H                   #DERIVS!C2DERTW!SOLVER I   #WRKC P   #WRK1 Q   #IY R                                                                                                                                         I     0                     #C2DERTW%SOLVER%GAIN J   #C2DERTW%SOLVER%UGAIN K   #C2DERTW%SOLVER%THETA L   #C2DERTW%SOLVER%ALPHA M   #C2DERTW%SOLVER%BETA N   #C2DERTW%SOLVER%DYDE O                                                    J             
                                                        K            
                                                        L            
                                                        M            
                                                        N             
                                                        O     (       
                                                           P      @              "    p A        p          p          p            p          p          p                                    D                                          Q      @              #    p A        p          p          p            p          p          p                                                                              R            #         @                                            S                   #DERIVS!CDERIV1D!SOLVER T   #F [   #DF \                                                                                                                                              T     0                     #CDERIV1D%SOLVER%GAIN U   #CDERIV1D%SOLVER%UGAIN V   #CDERIV1D%SOLVER%THETA W   #CDERIV1D%SOLVER%ALPHA X   #CDERIV1D%SOLVER%BETA Y   #CDERIV1D%SOLVER%DYDE Z                                                    U             
                                                        V            
                                                        W            
                                                        X            
                                                        Y             
                                                        Z     (       
                                                           [                   
 $    p          p            p                                    D                                          \                   
 %    p          p            p                                 ,      fn#fn    Ì   /      GRADV #   û        DERIVS!GRADV!WAVES !     ¬      GRADV%WAVES%WAVX !   >  ¬      GRADV%WAVES%WAVZ    ê  ¬      GRADV%WAVES%C      Ü   a   GRADV%U    r  Ü   a   GRADV%V    N  Ü   a   GRADV%W    *  Ü   a   GRADV%U11      Ü   a   GRADV%U12    â  Ü   a   GRADV%U13    ¾	  Ü   a   GRADV%U21    
  Ü   a   GRADV%U22    v  Ü   a   GRADV%U23    R  Ü   a   GRADV%U31    .  Ü   a   GRADV%U32    
  Ü   a   GRADV%U33    æ  Ü   a   GRADV%LU    Â  Ü   a   GRADV%LV      Ü   a   GRADV%LW    z  Ç       CDERIV %   A  ë      DERIVS!CDERIV!SOLVER #   ,  P      CDERIV%SOLVER%GAIN $   |  P      CDERIV%SOLVER%UGAIN $   Ì  P      CDERIV%SOLVER%THETA $     P      CDERIV%SOLVER%ALPHA #   l  P      CDERIV%SOLVER%BETA #   ¼  P      CDERIV%SOLVER%DYDE      Ü   a   CDERIV%F    è  Ü   a   CDERIV%DF    Ä  ¿       C0DERBW      Ü   a   C0DERBW%WRKC    _  Ü   a   C0DERBW%WRK1    ;  H   a   C0DERBW%IY      Ú       C1DERBW &   ]  ñ      DERIVS!C1DERBW!SOLVER $   N  P      C1DERBW%SOLVER%GAIN %     P      C1DERBW%SOLVER%UGAIN %   î  P      C1DERBW%SOLVER%THETA %   >  P      C1DERBW%SOLVER%ALPHA $     P      C1DERBW%SOLVER%BETA $   Þ  P      C1DERBW%SOLVER%DYDE    .  Ü   a   C1DERBW%WRKC    
  Ü   a   C1DERBW%WRK1    æ  H   a   C1DERBW%IY    .  Ú       C2DERBW &      ñ      DERIVS!C2DERBW!SOLVER $   ù   P      C2DERBW%SOLVER%GAIN %   I!  P      C2DERBW%SOLVER%UGAIN %   !  P      C2DERBW%SOLVER%THETA %   é!  P      C2DERBW%SOLVER%ALPHA $   9"  P      C2DERBW%SOLVER%BETA $   "  P      C2DERBW%SOLVER%DYDE    Ù"  Ü   a   C2DERBW%WRKC    µ#  Ü   a   C2DERBW%WRK1    $  H   a   C2DERBW%IY    Ù$  ¿       C0DERTW    %  Ü   a   C0DERTW%WRKC    t&  Ü   a   C0DERTW%WRK1    P'  H   a   C0DERTW%IY    '  Ú       C1DERTW &   r(  ñ      DERIVS!C1DERTW!SOLVER $   c)  P      C1DERTW%SOLVER%GAIN %   ³)  P      C1DERTW%SOLVER%UGAIN %   *  P      C1DERTW%SOLVER%THETA %   S*  P      C1DERTW%SOLVER%ALPHA $   £*  P      C1DERTW%SOLVER%BETA $   ó*  P      C1DERTW%SOLVER%DYDE    C+  Ü   a   C1DERTW%WRKC    ,  Ü   a   C1DERTW%WRK1    û,  H   a   C1DERTW%IY    C-  Ú       C2DERTW &   .  ñ      DERIVS!C2DERTW!SOLVER $   /  P      C2DERTW%SOLVER%GAIN %   ^/  P      C2DERTW%SOLVER%UGAIN %   ®/  P      C2DERTW%SOLVER%THETA %   þ/  P      C2DERTW%SOLVER%ALPHA $   N0  P      C2DERTW%SOLVER%BETA $   0  P      C2DERTW%SOLVER%DYDE    î0  Ü   a   C2DERTW%WRKC    Ê1  Ü   a   C2DERTW%WRK1    ¦2  H   a   C2DERTW%IY    î2  Ó       CDERIV1D '   Á3  ÷      DERIVS!CDERIV1D!SOLVER %   ¸4  P      CDERIV1D%SOLVER%GAIN &   5  P      CDERIV1D%SOLVER%UGAIN &   X5  P      CDERIV1D%SOLVER%THETA &   ¨5  P      CDERIV1D%SOLVER%ALPHA %   ø5  P      CDERIV1D%SOLVER%BETA %   H6  P      CDERIV1D%SOLVER%DYDE    6     a   CDERIV1D%F    47     a   CDERIV1D%DF 