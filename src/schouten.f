c     Schouten HJA. "Measuring pairwise interobserver agreement when all subjects are judged by the same observers." Statistica Neerlandica 1982; 36: 45-61.

c     Implementation by:
c     Mark Clements <mark.clements@ki.se>
c     GPL-3 or higher

c     Note: there was an error if I tried to remove some of the
c     subroutine parameters (e.g. pab..wa). Valgrind did not show any
c     evidence of a memory leak.

      subroutine schouten(N,M,L,data,w,kab,ka,kappa,pab,pa,p,ma,qab,qa
     $     ,q,oab,eab,oa,ea,o,e,wa,varkab,varka,vark,covkka,chi,var0kab
     $     ,var0ka,var0k,wab)
c     implicit none
      integer n, m, l, data(N,M)
      double precision w(L,L)
      double precision pab(M,M,L,L), pa(M,L,L), p(L,L)
      double precision ma(M,L), qab(M,M,L,L), qa(M,L,L)
      double precision q(L,L), oab(M,M), eab(M,M)
      double precision kab(M,M), oa(M), ea(M), ka(M)
      double precision o,e,kappa
      double precision wa(M,L),dab(M,M,L,L),sab(M,M)
      double precision varkab(M,M)
      double precision sa(M), varka(M), s, vark, covkka(M)
      double precision drow,darow,chi(M),wab(M,M)
      double precision var0kab(M,M),var0ka(M),var0k
      integer row,i,j,a,b,a2

      do 10 row=1,N
         do 20 a=1,M
            do 30 b=1,M
               pab(a,b,data(row,a),data(row,b)) = pab(a,b,data(row,a)
     $              ,data(row,b)) + 1.0/N
               
 30         continue
 20      continue
 10   continue

      do 110 a=1,M
         do 120 b=1,M
            do 130 i=1,L
               do 140 j=1,L
                  if(a .ne. b) then
                     pa(a,i,j) = pa(a,i,j)+pab(a,b,i,j)/(M-1)
                  end if
 140           continue
 130        continue
 120     continue
 110  continue

      do 210 a=1,M
         do 230 i=1,L
            do 240 j=1,L
               p(i,j) = p(i,j)+pa(a,i,j)/M
               ma(a,i) = ma(a,i)+pa(a,i,j)
 240        continue
 230     continue
 210  continue

      do 310 a=1,M
         do 320 b=1,M
            do 330 i=1,L
               do 340 j=1,L
                  qab(a,b,i,j) = ma(a,i)*ma(b,j)
                  if(a .ne. b) then
                     qa(a,i,j) = qa(a,i,j)+qab(a,b,i,j)/(M-1)
                  end if
 340           continue
 330        continue
 320     continue
 310  continue

      do 510 a=1,M
         do 530 i=1,L
            do 540 j=1,L
               q(i,j) = q(i,j)+qa(a,i,j)/M
 540        continue
 530     continue
 510  continue

      do 610 a=1,M
         do 620 b=1,M
            do 630 i=1,L
               do 640 j=1,L
                  oab(a,b) = oab(a,b)+pab(a,b,i,j)*w(i,j)
                  eab(a,b) = eab(a,b)+qab(a,b,i,j)*w(i,j)
 640           continue
 630        continue
            kab(a,b) = (oab(a,b)-eab(a,b))/(1.0-eab(a,b))
 620     continue
 610  continue

      do 710 a=1,M
         do 720 b=1,M
            if(a .ne. b) then
               oa(a) = oa(a)+oab(a,b)/(M-1)
               ea(a) = ea(a)+eab(a,b)/(M-1)
            end if
 720     continue
         ka(a) = (oa(a)-ea(a))/(1.0-ea(a))
         o = o+oa(a)/M
         e = e+ea(a)/M
 710  continue
      kappa = (o-e)/(1.0-e)

c     Variance calculations

c     Matrix multiplication
      call dgemm('n','n',M,L,L,1.0d0,ma,M,w,L,0.0d0,wa,M)

      do 910 a=1,M
         do 920 b=1,M
            do 930 i=1,L
               do 940 j=1,L
                  dab(a,b,i,j) = (1.0-eab(a,b))*w(i,j) - (1.0-oab(a,b))
     $                 *(wa(b,i)+wa(a,j))
                  varkab(a,b) = varkab(a,b)+pab(a,b,i,j)*dab(a,b,i,j)**2
     $                 /N/(1.0-eab(a,b))**4
 940           continue
 930        continue
            sab(a,b) = eab(a,b)*oab(a,b)-2*eab(a,b)+oab(a,b)
            varkab(a,b)=varkab(a,b)-sab(a,b)**2/N/(1.0-eab(a,b))**4
 920     continue
 910  continue

      do 1010 a=1,M
         sa(a) = ea(a)*oa(a)-2.0*ea(a)+oa(a)
         do 1020 row=1,N
            darow = 0.0
            do 1030 b=1,M
               if (a .ne. b) then
                  darow = darow + (1.0-ea(a))*w(data(row,a)
     $                 ,data(row,b))/(M-1) - (1.0-oa(a))*(wa(b,data(row
     $                 ,a))+wa(a,data(row,b)))/(M-1)
               end if
 1030       continue
            varka(a) = varka(a) + darow**2/N
 1020    continue
         varka(a) = (varka(a) - sa(a)**2) / N / (1.0-ea(a))**4
 1010 continue
            
      s = e*o-2.0*e+o
      do 1110 row=1,N
         drow = 0.0
         do 1120 a=1,M
            do 1130 b=1,M
               if (a .ne. b) then
                  drow = drow + (1.0-e)*w(data(row,a)
     $                 ,data(row,b))/M/(M-1) - 2.0*(1.0-o)*wa(b,data(row
     $                 ,a))/M/(M-1)
               end if
 1130       continue
 1120    continue
         vark = vark + drow**2/N
 1110 continue
      vark = (vark - s**2) / N / (1.0-e)**4

c     Covariance calculations      

      do 1210 a=1,M
         do 1220 row=1,N
            drow = 0.0
            darow = 0.0
            do 1240 b=1,M
               if (a .ne. b) then
                  darow = darow + (1.0-ea(a))*w(data(row,a)
     $                 ,data(row,b))/(M-1) - (1.0-oa(a))*(wa(b,data(row
     $                 ,a))+wa(a,data(row,b)))/(M-1)
               end if
c     Inefficient inner loop
               do 1250 a2=1,M
                  if (a2 .ne. b) then
                     drow = drow + (1.0-e)*w(data(row,a2) ,data(row,b))
     $                    /M/(M-1) - 2.0*(1.0-o)*wa(b,data(row ,a2))/M
     $                    /(M-1)
                  end if
 1250          continue
 1240       continue
            covkka(a) = covkka(a) + drow*darow/N
 1220    continue
         covkka(a) = (covkka(a) - sa(a)*s)/N/(1.0-e)**2/(1-ea(a))**2
         chi(a) = (kappa - ka(a))**2/(vark+varka(a)-2*covkka(a))
 1210 continue

c     Null hyptheses
      
      do 1310 a=1,M
         do 1320 b=1,M
            do 1330 i=1,L
               do 1340 j=1,L
                  wab(a,b) = wab(a,b) + qab(a,b,i,j)*(w(i,j)-wa(a,j)
     $                 -wa(b,i))**2
 1340          continue
 1330       continue
            wab(a,b) = wab(a,b) - eab(a,b)**2
            var0kab(a,b) = wab(a,b)/N/(1-eab(a,b))**2
            if (a .ne. b) then
               var0ka(a) = var0ka(a) + wab(a,b)/N/(M-1)**2/(1-ea(a))**2
               var0k = var0k + 2.0*wab(a,b)/N/M**2/(M-1)**2/(1-e)**2
            end if
 1320    continue
 1310 continue

      return
      end
