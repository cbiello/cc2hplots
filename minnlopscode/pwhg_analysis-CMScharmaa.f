
c     when using "dsig(1:rwl_num_weights)", i.e
c     rwl_weights(1:rwl_num_weights), we are actually going in real*4
c     precision. This is why some digits will differ wrt when we use
c     dsig0


      subroutine init_hist
      implicit none
      include 'pwhg_math.h'
      real * 8 dy,dpt
      real * 8 ymin, ymax, ptmin, ptmax
      integer icut
      integer, parameter :: ncuts=5
      real *8 jetcut(ncuts)
      character * 10 suffix(ncuts) 
      common/jcut/jetcut,suffix
      real *8 powheginput

      call inihists

cccc CB: Tiziano's settings
      dpt=5d0
      ptmin = 0d0
      ptmax = 200d0
      dy=0.4d0
      ymin = -3d0
      ymax = 3d0      
ccccccccccccccccccccccccccc

c ------------------------------------------------------
c NO CUTS
      call bookupeqbins('inc-integrated',1d0,0d0,1d0)
      call bookupeqbins('inc-ptH',dpt,ptmin,ptmax)
      call bookupeqbins('inc-ptc',dpt,ptmin,ptmax)
      call bookupeqbins('inc-yH',dy,ymin,ymax)
c ------------------------------------------------------
c aa cuts
      call bookupeqbins('fidaa-integrated',1d0,0d0,1d0)
      call bookupeqbins('fidaa-ptH',dpt,ptmin,ptmax)
      call bookupeqbins('fidaa-ptc',dpt,ptmin,ptmax)
      call bookupeqbins('fidaa-yH',dy,ymin,ymax)
c ------------------------------------------------------
c c-jet and aa cuts
      call bookupeqbins('fid-integrated',1d0,0d0,1d0)
      call bookupeqbins('fid-ptH',dpt,ptmin,ptmax)
      call bookupeqbins('fid-ptc',dpt,ptmin,ptmax)
      call bookupeqbins('fid-yH',dy,ymin,ymax)
c ------------------------------------------------------
c NO CUTS with BR included
      call bookupeqbins('incBR-integrated',1d0,0d0,1d0)
      call bookupeqbins('incBR-ptH',dpt,ptmin,ptmax)
      call bookupeqbins('incBR-ptc',dpt,ptmin,ptmax)
      call bookupeqbins('incBR-yH',dy,ymin,ymax)
c ------------------------------------------------------
c c-jet and aa cuts with BR included
      call bookupeqbins('fidaaBR-integrated',1d0,0d0,1d0)
      call bookupeqbins('fidaaBR-ptH',dpt,ptmin,ptmax)
      call bookupeqbins('fidaaBR-ptc',dpt,ptmin,ptmax)
      call bookupeqbins('fidaaBR-yH',dy,ymin,ymax)
c ------------------------------------------------------
c c-jet and aa cuts with BR included
      call bookupeqbins('fidBR-integrated',1d0,0d0,1d0)
      call bookupeqbins('fidBR-ptH',dpt,ptmin,ptmax)
      call bookupeqbins('fidBR-ptc',dpt,ptmin,ptmax)
      call bookupeqbins('fidBR-yH',dy,ymin,ymax)

ccccc TODO: Haa only CUTS
      
      end

      
     
      subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'pwhg_weights.h'
      include 'nlegborn.h'
      include 'pwhg_rad.h'
c      include 'pwhg_lhrwgt.h'
      include 'pwhg_rwl.h'
      logical ini
      data ini/.true./
      save ini
C     allow multiweights 
      real * 8 dsig0,dsig(1:weights_max)
c     tell to the analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer, parameter :: ncuts=5
      real *8 jetcut(ncuts)
      character * 10 suffix(ncuts) 
      common/jcut/jetcut,suffix
C     Higgs variables 
      real*8 pH(4),mH,ptH,yH,etaH
c     arrays to reconstruct jets
      integer maxjet
      parameter (maxjet=2048)
      real *8 ptmin
      real *8  ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1     phij(maxjet),pj(4,maxjet),rr,ptrel(4),
     2     yj1,etaj1,ptj1,mj1,dyhj1,detahj1,dphihj1,drhj1
      integer j1,found,mjets,ihep,icut,i
      real *8 powheginput
      real *8 corrfactor
      real *8 ptHcut, ptJcut, ptHmax, ptHmin
      real *8 ptsave
      real *8 dyhj2,detahj2,dphihj2,drhj2
      integer ijet, ihardest, isecond
      logical condition
      real *8 BRHaa,mgamgam
      real *8 dya1a2,detaa1a2,dphia1a2,dra1a2
      real *8 eta_photon1,eta_photon2,eta_photon_a,eta_photon_b
      real *8 m_photon1,m_photon2,m_photon_a,m_photon_b
      real *8 y_photon1,y_photon2,y_photon_a,y_photon_b
      real *8 pt_photon1,pt_photon2,pt_photon_a,pt_photon_b
      real *8 mHiggs
      parameter (mHiggs=125.0d0)
      integer mu, n_photon
      real *8 ptgamgam
      integer photonvec(3)
      real *8 dsigBR(1:weights_max)
      integer total_c_jets,counter_c_jets,counter_all_jets
      real *8 kt_c_jet(maxjet),eta_c_jet(maxjet),rap_c_jet(maxjet), phi_c_jet(maxjet)
      real *8 pt_cjet_cut_CMS, eta_cjet_cut_CMS,pt_a1_cut_CMS,pt_a2_cut_CMS, eta_a_cut_CMS
      integer n_cjets
      real *8 p_c_jets(4,maxjet), p_all_jets(4,maxjet)
      real *8 etamax, tmp
      integer j
      logical is_cjet_array(maxjet)


      BRHaa=0.00227d0 !Branching fraction of H->gam+gam decay
      
      corrfactor=1d0
      if(powheginput("#btildeviol").eq.1.and.WHCPRG.ne.'NLO') then
         if(rad_type.eq.1) then
            corrfactor=powheginput('corr_btilde')
         elseif(rad_type.eq.2) then
            corrfactor=powheginput('corr_remnant')
         else
            print*, 'no rew'
         endif
      endif


      if (ini) then
         write(*,*) '*****************************'
         if(whcprg.eq.'NLO') then
            write(*,*) '       NLO ANALYSIS'
            weights_num=0
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE ANALYSIS'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         elseif(WHCPRG.eq.'PY8   ') then
            write (*,*) '           PYTHIA 8 ANALYSIS        '
         endif
         write(*,*) '*****************************'
                  
         write(*,*) ''
         write(*,*) '*****************************'
         write(*,*) '** weights_num     = ',weights_num
         write(*,*) '** rwl_num_weights = ',rwl_num_weights
         write(*,*) '** rwl_num_groups = ',rwl_num_groups
         write(*,*) '*****************************'
         write(*,*) ''
            
         if(weights_num.eq.0.and.rwl_num_weights.eq.0) then
            call setupmulti(1)
         else if(weights_num.ne.0.and.rwl_num_weights.eq.0) then
            call setupmulti(weights_num)
         else if(weights_num.eq.0.and.rwl_num_weights.ne.0) then
            call setupmulti(rwl_num_weights)
         else
            call setupmulti(rwl_num_weights)
         endif
            
         if(weights_num.eq.0.and.rwl_num_weights.gt.weights_max) then
            write(*,*) 'ERROR:'
            write(*,*) 'incoming number of weights (rwl_num_weights)'
            write(*,*) 'is greater than declared dsig and bWdsig    '
            write(*,*) 'array length.'
            stop
         endif

         ini=.false.
      endif

      dsig=0
      
      if(weights_num.eq.0.and.rwl_num_weights.eq.0) then
         dsig(1)=dsig0
      else if(weights_num.ne.0.and.rwl_num_weights.eq.0) then
         dsig(1:weights_num)=weights_val(1:weights_num)
      else if(weights_num.eq.0.and.rwl_num_weights.ne.0) then
         dsig(1:rwl_num_weights)=rwl_weights(1:rwl_num_weights)
      else
         dsig(1:rwl_num_weights)=rwl_weights(1:rwl_num_weights)
      endif
      
      if(sum(abs(dsig)).eq.0) return

      dsig=dsig*corrfactor

      do i=1,rwl_num_weights
        if(abs(dsig(i))>1d2 .or. dsig(i)+1 .eq. dsig(i)) then
          write(*,*) "LARGE weight. DISCARDING EVENT, i, weight = ",i, dsig(i)
          return
        endif
      enddo
      
c     Loop over final state particles to find Higgs
      found=0
      if(WHCPRG.eq.'PY8   ') then
         n_photon=0
         do ihep=1,nhep
            if((idhep(ihep).eq.25).and.(isthep(ihep).eq.1)) then !'stable' Higgs
                  pH(1:4) = phep(1:4,ihep)
                  found=1
            endif  
           if(idhep(ihep).eq.22) then
            if((isthep(ihep).eq.1).and.(idhep(jmohep(1,ihep)).eq.25)) then
               n_photon=n_photon+1
               photonvec(n_photon)=ihep
            endif
            endif
         enddo !nhep
      else
         do ihep=1,nhep
            if (((isthep(ihep).eq.1).or.(isthep(ihep).eq.2)
     #.or.(isthep(ihep).eq.155).or.(isthep(ihep).eq.195))
     #.and.(idhep(ihep).eq.25)) then
               j1=ihep
               found=found+1
               pH(1:4) = phep(1:4,j1)
            endif
         enddo
      endif
      
c     Check if there is still stable Higgs after PY8 showering
      if (WHCPRG.eq.'PY8') then                                                                                                                 
         if(found.gt.0) then
            write(*,*) 'CB: Stable Higgs in the event.'
            write(*,*) 'Why the Higgs did not decay into two photons?'
            call exit(-1)
         else
            dsigBR=dsig*BRHaa
         endif
      endif
      if(n_photon.ne.2) then
            write(*,*) 'ERROR: Higgs not found and less than 2 photons'
            call exit(1)
      else
            call getyetaptmass(phep(:,photonvec(1)),y_photon_a,eta_photon_a,pt_photon_a,m_photon_a)
            call getyetaptmass(phep(:,photonvec(2)),y_photon_b,eta_photon_b,pt_photon_b,m_photon_b)
            call getdydetadphidr(phep(:,photonvec(1)),phep(:,photonvec(2)),dya1a2,detaa1a2,dphia1a2,dra1a2)

            if(pt_photon_a .gt. pt_photon_b) then
               y_photon1   = y_photon_a
               eta_photon1 = eta_photon_a
               pt_photon1  = pt_photon_a
               m_photon1   = m_photon_a

               y_photon2   = y_photon_b
               eta_photon2 = eta_photon_b
               pt_photon2  = pt_photon_b
               m_photon2   = m_photon_b
            else
               y_photon1   = y_photon_b
               eta_photon1 = eta_photon_b
               pt_photon1  = pt_photon_b
               m_photon1   = m_photon_b

               y_photon2   = y_photon_a
               eta_photon2 = eta_photon_a
               pt_photon2  = pt_photon_a
               m_photon2   = m_photon_a
            endif
c  reconstruct Higgs momentum from hardest and second hardest photon (only two photons should be there)
            do mu=1,4
               pH(mu)=phep(mu,photonvec(1))+phep(mu,photonvec(2))
            enddo
            mgamgam=dsqrt(abs(pH(4)**2-pH(1)**2-pH(2)**2-pH(3)**2))
            ptgamgam=dsqrt(abs(pH(1)**2+pH(2)**2))
            call getyetaptmass(pH,yH,etaH,ptH,mH)
      endif

      call filld('inc-integrated', 0.5d0, dsig)
      call filld('inc-ptH', ptH, dsig)
      call filld('inc-yH', yH, dsig)
      call filld('incBR-integrated', 0.5d0, dsigBR)
      call filld('incBR-ptH', ptH, dsigBR)
      call filld('incBR-yH', yH, dsigBR)
      
cccccccccccccccc-JET ANALYSIS cccccccccccccc
      rr=0.4d0
      ptmin=0d0
      etamax=10d0
c      call buildjets(1,rr,ptmin,mjets,ktj,etaj,rapj,phij,ptrel,pj)
      call buildjetsandcjets(1,rr,ptmin,mjets,ktj,etaj,rapj,phij,ptrel,pj,is_cjet_array,total_c_jets)

      if(mjets.eq.0) then
         print*, 'Error in pwhg_analysis! There must be at least 1j'
         stop
      endif

      counter_c_jets=0
      counter_all_jets=0

      do i=1,mjets
         if(dabs(etaj(i)).lt.etamax) then
            if(is_cjet_array(i)) then
               counter_c_jets=counter_c_jets+1
               p_c_jets(:,counter_c_jets) = pj(:,i)
            endif
            counter_all_jets=counter_all_jets+1
            p_all_jets(:,counter_all_jets) = pj(:,i)
         endif
      enddo
      do j=1,counter_c_jets
         call getyetaptmass2(p_c_jets(:,j),rap_c_jet(j),eta_c_jet(j),kt_c_jet(j),tmp)
         phi_c_jet(j)=atan2(p_c_jets(2,j),p_c_jets(1,j))
      enddo

      ptmin=0d0
      etamax=10d0
      ihardest=-1
      if(counter_c_jets.ge.1) then
c     if at least a jet was found, loop over jets and find the hardest within rapidity cut
         ptsave=-1d0
         do ijet=1,counter_c_jets
            condition=(dabs(eta_c_jet(ijet)).le.etamax).and.(kt_c_jet(ijet).ge.ptsave)
            if(condition) then
               ihardest=ijet
               ptsave=kt_c_jet(ihardest)
            endif
            if(kt_c_jet(ijet).lt.ptmin) then
               write(*,*) 'ERROR1: this cannot happen'
            endif
         enddo
c     if a good jet is found, fill histograms
         if(ihardest.ne.-1) then
            call filld('inc-ptc',kt_c_jet(ihardest),dsig)
            call filld('incBR-ptc',kt_c_jet(ihardest),dsigBR)
         endif
      endif
      
c CMS CUTS
      pt_cjet_cut_CMS = 25d0
      eta_cjet_cut_CMS = 2.5d0
      pt_a1_cut_CMS = 125d0*0.33d0
      pt_a2_cut_CMS = 125d0*0.25d0
      eta_a_cut_CMS = 2.5d0
c Search for aa satisfying the cuts
      if(abs(eta_photon1) .lt.eta_a_cut_CMS .and. pt_photon1 .gt. pt_a1_cut_CMS) then
            if(abs(eta_photon2) .lt.eta_a_cut_CMS .and. pt_photon2 .gt. pt_a2_cut_CMS) then

                call filld('fidaa-integrated', 0.5d0, dsig)
                call filld('fidaa-ptH', ptH, dsig)
                call filld('fidaa-ptc', kt_c_jet(ihardest), dsig)
                call filld('fidaa-yH', yH, dsig)
                call filld('fidaaBR-integrated', 0.5d0, dsigBR)
                call filld('fidaaBR-ptH', ptH, dsigBR)
                call filld('fidaaBR-ptc', kt_c_jet(ihardest), dsigBR)
                call filld('fidaaBR-yH', yH, dsigBR)

            endif
      endif
c Search for cjet satisfying the cuts
      if(abs(eta_c_jet(ihardest)).lt.eta_cjet_cut_CMS .and. kt_c_jet(ihardest).gt. pt_cjet_cut_CMS) then
         if(abs(eta_photon1) .lt.eta_a_cut_CMS .and. pt_photon1 .gt. pt_a1_cut_CMS) then
            if(abs(eta_photon2) .lt.eta_a_cut_CMS .and. pt_photon2 .gt. pt_a2_cut_CMS) then

                call filld('fid-integrated', 0.5d0, dsig)
                call filld('fid-ptH', ptH, dsig)
                call filld('fid-ptc', kt_c_jet(ihardest), dsig)
                call filld('fid-yH', yH, dsig)
                call filld('fidBR-integrated', 0.5d0, dsigBR)
                call filld('fidBR-ptH', ptH, dsigBR)
                call filld('fidBR-ptc', kt_c_jet(ihardest), dsigBR)
                call filld('fidBR-yH', yH, dsigBR)

            endif
         endif
      endif
 
      end
      


      subroutine buildjets(iflag,rr,ptmin,mjets,kt,eta,rap,phi,
     $     ptrel,pjet)
c     arrays to reconstruct jets, radius parameter rr
      implicit none
c     tell to the analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer iflag,mjets
      real * 8  rr,ptmin,kt(*),eta(*),rap(*),
     1     phi(*),ptrel(3),pjet(4,*)
      include   'hepevt.h'
      include  'LesHouches.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu,i
      real * 8 r,palg,tmp
      logical islept
      external islept
      real * 8 vec(3),pjetin(0:3),pjetout(0:3),beta,
     $     ptrackin(0:3),ptrackout(0:3)
      real * 8 get_ptrel
      external get_ptrel
C - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
      enddo      
      ntracks=0
      do j=1,maxjet
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
      enddo
      if(iflag.eq.1) then
C     - Extract final state particles to feed to jet finder
         if(WHCPRG.eq.'PY8   ') then
            do j=1,nhep
               if((isthep(j).eq.1).and..not.idhep(j).eq.25) then
                  if(ntracks.eq.maxtrack) then
                     write(*,*) 'analyze: need to increase maxtrack!'
                     write(*,*) 'ntracks: ',ntracks
                     stop
                  endif
                  ntracks=ntracks+1
                  do mu=1,4
                     ptrack(mu,ntracks)=phep(mu,j)
                  enddo
c                  write(13,*) isthep(j),ptrack(:,ntracks) 
                  itrackhep(ntracks)=j
               else
c                  write(14,*) isthep(j),ptrack(:,ntracks) 
               endif
            enddo
c            stop
cc Stop fr debugging
         else
            do j=1,nhep
c     all but the Higgs
               if (isthep(j).eq.1.and..not.idhep(j).eq.25) then
                  if(ntracks.eq.maxtrack) then
                     write(*,*) 'analyze: need to increase maxtrack!'
                     write(*,*) 'ntracks: ',ntracks
                     stop
                  endif
                  ntracks=ntracks+1
                  do mu=1,4
                     ptrack(mu,ntracks)=phep(mu,j)
                  enddo
                  itrackhep(ntracks)=j
               endif
            enddo
         endif
      else
         do j=1,nup
            if (istup(j).eq.1.and..not.islept(idup(j))) then
               if(ntracks.eq.maxtrack) then
                  write(*,*) 'analyze: need to increase maxtrack!'
                  write(*,*) 'ntracks: ',ntracks
                  stop
               endif
               ntracks=ntracks+1
               do mu=1,4
                  ptrack(mu,ntracks)=pup(mu,j)
               enddo
               itrackhep(ntracks)=j
            endif
         enddo
      endif
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C --------------------------------------------------------------------- C
C     R = 0.7   radius parameter
c palg=1 is standard kt, -1 is antikt
      palg=-1
      r=rr
c      ptmin=20d0 
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec)
      mjets=njets
      if(njets.eq.0) return
c check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4) then
         write(*,*) ' bug!'
      endif
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         call getyetaptmass(pjet(:,j),rap(j),eta(j),kt(j),tmp)
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo

c     loop over the hardest 3 jets
      do j=1,min(njets,3)
         do mu=1,3
            pjetin(mu) = pjet(mu,j)
         enddo
         pjetin(0) = pjet(4,j)         
         vec(1)=0d0
         vec(2)=0d0
         vec(3)=1d0
         beta = -pjet(3,j)/pjet(4,j)
         call mboost(1,vec,beta,pjetin,pjetout)         
c     write(*,*) pjetout
         ptrel(j) = 0
         do i=1,ntracks
            if (jetvec(i).eq.j) then
               do mu=1,3
                  ptrackin(mu) = ptrack(mu,i)
               enddo
               ptrackin(0) = ptrack(4,i)
               call mboost(1,vec,beta,ptrackin,ptrackout) 
               ptrel(j) = ptrel(j) + get_ptrel(ptrackout,pjetout)
            endif
         enddo
      enddo
      end

C     a number of handy functions used by generic analyses 

      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real *8 tiny
      parameter (tiny=1.d-5)
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt.lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      else
         eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      mass=sqrt(abs(p(4)**2-pv**2))
      end

      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(*),p2(*),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call getyetaptmass(p1,y1,eta1,pt1,mass1)
      call getyetaptmass(p2,y2,eta2,pt2,mass2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end

      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.15) then
         islept = .true.
      else
         islept = .false.
      endif
      end

      function get_ptrel(pin,pjet)
      implicit none
      real * 8 get_ptrel,pin(0:3),pjet(0:3)
      real * 8 pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      get_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end


c     ROUTINE FOR C-tagging. I call the c-jets as b-jet in this routine
c     beacuse copy and paste is easy
      subroutine buildjetsandcjets(iflag,rr,ptmin,mjets,kt,eta,rap,phi,
     $     ptrel,pjet,is_bjet_array,total_b_jets)
c     arrays to reconstruct jets, radius parameter rr
      implicit none
c     tell to the analysis file which program is running it
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer iflag,mjets
      real * 8  rr,ptmin,kt(*),eta(*),rap(*),
     1     phi(*),ptrel(3),pjet(4,*)
      logical is_bjet_array(*)
      include   'hepevt.h'
      include  'LesHouches.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu,i
      real * 8 r,palg,tmp
      real * 8 vec(3),pjetin(0:3),pjetout(0:3),beta,
     $     ptrackin(0:3),ptrackout(0:3)
      real * 8 get_ptrel
      external get_ptrel
      logical is_C_hadron,is_CBAR_hadron,is_neutrino
      external is_C_hadron,is_CBAR_hadron,is_neutrino
      logical is_b_track(maxtrack)
      integer const_indices(maxtrack),nconst,ijet
      integer nbjet_array(maxjet),
     $     nbbarjet_array(maxjet),jetinfo(maxjet),id,nb,
     $     nbbar,nbjet,nbbarjet,total_b_jets
C - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
      enddo      
      ntracks=0
      do j=1,maxjet
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
         is_bjet_array(j) = .false. 
      enddo
      if(iflag.eq.1) then
C     - Extract final state particles to feed to jet finder
         if(WHCPRG.eq.'PY8   ') then
            do j=1,nhep
c all but top/anti-top
c     exclude leptons, gauge and higgs bosons, but include gluons             
               if(isthep(j).eq.1.and.((abs(idhep(j)).le.5.or.abs(idhep(j)).ge.40
               ! JM: changed isthep == 1 by is_final_state, check! (see also further changes below)
               !if(is_final_state(j).and.((abs(idhep(j)).le.4.or.abs(idhep(j)).ge.40
     &              .or.abs(idhep(j)).eq.21))) then
                  if(ntracks.eq.maxtrack) then
                     write(*,*) 'analyze: need to increase maxtrack!'
                     write(*,*) 'ntracks: ',ntracks
                     stop
                  endif
                  ntracks=ntracks+1
                  do mu=1,4
                     ptrack(mu,ntracks)=phep(mu,j)
                  enddo
                  itrackhep(ntracks)=j
               endif
            enddo
         else
           do j=1,nhep
c all but top/anti-top
               if(isthep(j).eq.1.and.((abs(idhep(j)).le.5.or.abs(idhep(j)).ge.40
     &              .or.abs(idhep(j)).eq.21))) then
                  if(ntracks.eq.maxtrack) then
                     write(*,*) 'analyze: need to increase maxtrack!'
                     write(*,*) 'ntracks: ',ntracks
                     stop
                  endif
                  ntracks=ntracks+1
                  do mu=1,4
                     ptrack(mu,ntracks)=phep(mu,j)
                  enddo
                  itrackhep(ntracks)=j
               endif
            enddo
         endif
      else
         do j=1,nup
            if(istup(j).eq.1.and.((abs(idhep(j)).le.5.or.abs(idhep(j)).ge.40
     &           .or.abs(idhep(j)).eq.21))) then
               if(ntracks.eq.maxtrack) then
                  write(*,*) 'analyze: need to increase maxtrack!'
                  write(*,*) 'ntracks: ',ntracks
                  stop
               endif
               ntracks=ntracks+1
               do mu=1,4
                  ptrack(mu,ntracks)=pup(mu,j)
               enddo
               itrackhep(ntracks)=j
            endif
         enddo
      endif
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C --------------------------------------------------------------------- C
C     R = r  radius parameter
c palg=1 is standard kt, -1 is antikt
      palg=-1
      r=rr
c      ptmin=20d0 
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $                        jetvec)
      mjets=njets

c JM: added      
c----------------------------------------------------------------------
c     find in which ptrack the B hadrons ended up
      nbjet_array = 0
      nbbarjet_array = 0
      nbjet=0
      nbbarjet=0
c     loop over tracks
      do i=1,ntracks
         id=idhep(itrackhep(i))
         if (is_C_hadron(id)) then
            nbjet=nbjet+1
            nbjet_array(nbjet)=jetvec(i)            
         elseif (is_CBAR_hadron(id)) then   
            nbbarjet=nbbarjet+1
            nbbarjet_array(nbbarjet)=jetvec(i)                        
         endif
      enddo
      
      total_b_jets = 0
      do i=1,njets
         jetinfo(i)=0
         is_bjet_array(i)=.false.
         do j=1,nbjet
            if (i.eq.nbjet_array(j)) then
               is_bjet_array(i)=.true.
            endif
         enddo
         do j=1,nbbarjet
            if (i.eq.nbbarjet_array(j)) then
               is_bjet_array(i)=.true.
             endif
          enddo
          if(is_bjet_array(i)) total_b_jets = total_b_jets + 1 
       enddo

c----------------------------------------------------------------------

      
      if(njets.eq.0) return
c check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4) then
         write(*,*) ' bug!'
      endif
C --------------------------------------------------------------------- C
C - Computing arrays of useful kinematics quantities for hardest jets - C
C --------------------------------------------------------------------- C
      do j=1,mjets
         call getyetaptmass2(pjet(:,j),rap(j),eta(j),kt(j),tmp)
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo

c     loop over the hardest 3 jets
      do j=1,min(njets,3)
         do mu=1,3
            pjetin(mu) = pjet(mu,j)
         enddo
         pjetin(0) = pjet(4,j)         
         vec(1)=0d0
         vec(2)=0d0
         vec(3)=1d0
         beta = -pjet(3,j)/pjet(4,j)
         call mboost(1,vec,beta,pjetin,pjetout)         
c     write(*,*) pjetout
         ptrel(j) = 0
         do i=1,ntracks
            if (jetvec(i).eq.j) then
               do mu=1,3
                  ptrackin(mu) = ptrack(mu,i)
               enddo
               ptrackin(0) = ptrack(4,i)
               call mboost(1,vec,beta,ptrackin,ptrackout) 
               ptrel(j) = ptrel(j) + get_ptrel(ptrackout,pjetout)
            endif
         enddo
      enddo
      end


      function is_C_hadron(id)
      implicit none
      logical is_C_hadron
      integer id
      is_C_hadron=((id.gt.400).and.(id.lt.500)).or.
     $     ((id.gt.4000).and.(id.lt.5000)).or.(id.eq.4)
      end

      function is_CBAR_hadron(id)
      implicit none
      logical is_CBAR_hadron
      integer id
      is_CBAR_hadron=((id.gt.-500).and.(id.lt.-400)).or.
     $     ((id.gt.-5000).and.(id.lt.-4000)).or.(id.eq.-4)
      end

      subroutine getyetaptmass2(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real *8 tiny
      parameter (tiny=1.d-5)
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt.lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      else
         eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      mass=sqrt(abs(p(4)**2-pv**2))
      end

