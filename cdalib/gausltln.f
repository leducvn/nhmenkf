C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
      SUBROUTINE GAUSLTLN (ALAT,ALON, IM,JM)                                    
                                                                                
*---------------------------------------------------------------------          
*  GAUSS格子の緯度,経度を計算する                                       
*                                  元祖  'A0568.NEW.FORT(GAUSS)'            
*                           CDA用に修正   1994.9.26   K.ONOGI               
*---------------------------------------------------------------------          
*  A ; COSINE OF COLATITUDE                                                     
*  W ; GAUSSIAN WEIGHT                                                          
*  JM; ORDER OF LEGENDRE FUNCTIONS                                              
*---------------------------------------------------------------------          
                                                                                
      PARAMETER (JMX = 320)                                                     
                                                                                
      IMPLICIT REAL*8 (A-H, O-Z)                                                
      REAL  *8 A(JMX), W(JMX)                                                   
      REAL  *4 ALAT(IM,JM), ALON(IM,JM)                                         
                                                                                
      LOGICAL LX
      DATA LX/.TRUE./                                                        
      IF (LX) THEN                                                              
*       配列の大きさのチェック                                              
        IF (JMX.LT.JM) THEN                                                     
          WRITE(6,*) 'GAUSLTLN 配列サイズ不適合'                            
          WRITE(6,*) 'JMX,JM=',JMX,JM                                           
          STOP 200                                                              
        END IF                                                                  
                                                                                
        ESP  = 1.D-14                                                           
        C    = (1.D0-(2.D0/3.14159265358979D0)**2)*0.25D0                       
        DRAD = 180.D0/3.1415926536D0                                            
                                                                                
        LX = .FALSE.                                                            
      END IF                                                                    
                                                                                
      FK  = JM                                                                  
      KK  = JM/2                                                                
      CALL BSSLZ1 (A,KK)                                                        
                                                                                
      DO 30 IS=1,KK                                                             
        ITER = 0                                                                
        XZ = COS(A(IS)/SQRT((FK+0.5D0)**2+C))                                   
   10   CONTINUE                                                                
        PKM2 = 1.D0                                                             
        PKM1 = XZ                                                               
        ITER = ITER + 1                                                         
        IF (ITER.GT.10) GOTO 700                                                
                                                                                
        DO 20 N=2,JM                                                            
          FN   = N                                                              
          PK   = ((2.D0*FN-1.D0)*XZ*PKM1-(FN-1.D0)*PKM2)/FN                     
          PKM2 = PKM1                                                           
          PKM1 = PK                                                             
   20   CONTINUE                                                                
                                                                                
        PKM1  = PKM2                                                            
        PKMRK = (FK*(PKM1-XZ*PK))/(1.D0-XZ**2)                                  
        SP    = PK/PKMRK                                                        
        XZ    = XZ - SP                                                         
        AVSP  = ABS(SP)                                                         
        IF (AVSP.GT.ESP) GOTO 10                                                
        A(IS) = XZ                                                              
        W(IS) = (2.D0*(1.D0-XZ**2))/(FK*PKM1)**2                                
   30 CONTINUE                                                                  
                                                                                
      IF( JM.EQ.KK*2 ) GOTO 50                                                  
      A(KK+1) = 0.D0                                                            
      PK = 2.D0/FK**2                                                           
                                                                                
      DO 40 N=2,JM,2                                                            
        FN = N                                                                  
        PK = PK*FN**2/(FN-1.D0)**2                                              
   40 CONTINUE                                                                  
      W(KK+1) = PK                                                              
                                                                                
   50 CONTINUE                                                                  
      DO 60 N=1,KK                                                              
        L    = JM+1-N                                                           
        A(L) = -A(N)                                                            
        W(L) =  W(N)                                                            
   60 CONTINUE                                                                  
                                                                                
      DO 70 I=1,IM                                                              
      DO 80 J=1,JM                                                              
        LL = JM+1-J                                                             
        ALAT(I,J) = DACOS(A(LL))*DRAD-90.D0                                     
        ALON(I,J) = 360.D0/FLOAT(IM) * (I-1)                                    
   80 CONTINUE                                                                  
   70 CONTINUE                                                                  
                                                                                
      RETURN                                                                    
                                                                                
  700 WRITE(6,6000)                                                             
 6000 FORMAT(//5X,14HERROR IN GAUAW//)                                          
                                                                                
      STOP                                                                      
      END                                                                       
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
      SUBROUTINE BSSLZ1 (BES, N)                                                
                                                                                
      IMPLICIT REAL*8(A-H,O-Z)                                                  
      REAL  *8 BES(N), BZ(50)                                                   
                                                                                
      DATA PI/3.14159265358979D0/                                               
      DATA BZ         / 2.4048255577D0, 5.5200781103D0,                         
     \  8.6537279129D0,11.7915344391D0,14.9309177086D0,18.0710639679D0,         
     \ 21.2116366299D0,24.3524715308D0,27.4934791320D0,30.6346064684D0,         
     \ 33.7758202136D0,36.9170983537D0,40.0584257646D0,43.1997917132D0,         
     \ 46.3411883717D0,49.4826098974D0,52.6240518411D0,55.7655107550D0,         
     \ 58.9069839261D0,62.0484691902D0,65.1899648002D0,68.3314693299D0,         
     \ 71.4729816036D0,74.6145006437D0,77.7560256304D0,80.8975558711D0,         
     \ 84.0390907769D0,87.1806298436D0,90.3221726372D0,93.4637187819D0,         
     \ 96.6052679510D0,99.7468198587D0,102.888374254D0,106.029930916D0,         
     \ 109.171489649D0,112.313050280D0,115.454612653D0,118.596176630D0,         
     \ 121.737742088D0,124.879308913D0,128.020877005D0,131.162446275D0,         
     \ 134.304016638D0,137.445588020D0,140.587160352D0,143.728733573D0,         
     \ 146.870307625D0,150.011882457D0,153.153458019D0,156.295034268D0/         
                                                                                
      NN = N                                                                    
      IF (N.LE.50) GOTO 12                                                      
      BES(50) = BZ(50)                                                          
                                                                                
      DO 5 J=51,N                                                               
        BES(J) = BES(J-1) + PI                                                  
    5 CONTINUE                                                                  
                                                                                
      NN=49                                                                     
                                                                                
   12 CONTINUE                                                                  
      DO 15 J=1,NN                                                              
        BES(J) = BZ(J)                                                          
   15 CONTINUE                                                                  
                                                                                
      RETURN                                                                    
      END                                                                       
