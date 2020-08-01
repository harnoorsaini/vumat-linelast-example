C     CONTENTS:
C           - Isotropic linear elastic material. Based on the example found 
C           in the Abaqus User Subroutines Reference Guide section 1.2.20
C           - The subroutine was only tested for 3D solid elements 
C     AUTHOR:
C           Harnoor Saini
C           Institute of Modelling and Simulation of Biomechanical Systems
C           University of Stuttgart
C           Stuttgart, Germany
C     CREATED:
C           11-10-2019
C     TESTED ON:
C           - Windows 10 64-bit v1803 OS build 17134.885
C           - Intel Core i7-6700HQ
C           - Abaqus 3DEXPERIENCE R2017x
C           - Intel Visual Fortran 64 Compiler V16.0.4.246 Build 20160811
C
C     NOTE:
C           - Abaqus generally does NOT use the "implicit none" statement, 
C           - therefore any variable starting with i, j, k, l, m and n are 
C           automatically treated as integers
C            
C     LICENSE:    
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C

      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock), cd(nblock), sti(nblock)
C
      character*80 cmname
C
      parameter( zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0,
     1  third = one/three, half = .5d0, twoThirds = two/three,
     2  threeHalfs = 1.5d0 )
C     
      integer i
      real e, xnu, twomu, sixmu, alamda, trace, stressPower, cd, sti
C
C     material parameters
      e     = props(1)
      xnu   = props(2)
      twomu  = e / ( one + xnu )
      sixmu  = three * twomu
      alamda = e * xnu / ( ( one + xnu ) * (one - two * xnu) )
C      
      do i = 1,nblock
C           update the stresses
            trace  = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
            stressNew(i,1) = stressOld(i,1) + alamda*trace + 
     1            twomu*strainInc(i,1)
            stressNew(i,2) = stressOld(i,2) + alamda*trace + 
     1            twomu*strainInc(i,2)
            stressNew(i,3) = stressOld(i,3) + alamda*trace + 
     1            twomu*strainInc(i,3)
            stressNew(i,4) = stressOld(i,4) + twomu*strainInc(i,4)
            stressNew(i,5) = stressOld(i,5) + twomu*strainInc(i,5)
            stressNew(i,6) = stressOld(i,6) + twomu*strainInc(i,6)
C           update the strain energy density
            stressPower = half * (
     1            (stressNew(i,1) * strainInc(i,1)) 
     2            + (stressNew(i,2) * strainInc(i,2))
     3            + (stressNew(i,3) * strainInc(i,3)) 
     4            + two * (stressNew(i,4) * strainInc(i,4)) 
     5            + two * (stressNew(i,5) * strainInc(i,5)) 
     6            + two * (stressNew(i,6) * strainInc(i,6)) )
            enerInternNew(i) = enerInternOld(i)
     1            + stressPower / density(i)
C           stable time increment estimation
            cd(i) = sqrt ( alamda + twomu / density(i) )
            sti(i) = charLength(i) / cd(i)
C           state variables
C           total strain in the 22 direction
            stateNew(i,1) = stateOld(i,1)+strainInc(i,2)
            stateNew(i,2) = enerInternNew(i)
            stateNew(i,3) = sti(i)
      end do
C
      return
      end