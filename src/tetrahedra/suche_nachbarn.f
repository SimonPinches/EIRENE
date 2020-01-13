

      SUBROUTINE EIRENE_SUCHE_NACHBARN

      USE EIRMOD_PARMMOD
      USE EIRMOD_CTETRA

      IMPLICIT NONE

      TYPE(TET_ELEM), POINTER :: CUR
C      TYPE(TET_ELEM), POINTER :: CUR2
      INTEGER :: ITET,IS,JTET,JS,
     .           IP1, i, j
C      INTEGER :: IC
      INTEGER :: JP(3),ip(3)
      INTEGER :: ITSIDE(3,4)
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/

      DO ITET=1,NTET      ! FOR ALL TETRAHEDRA
        DO IS=1,4         ! AND FOR ALL SIDES OF EACH TETRAHEDRON
          IF (NTBAR(IS,ITET) == 0) THEN   ! IF IT HAS NO NEIGHBOR YET
            IP(1)=NTECK(ITSIDE(1,IS),ITET)
            IP(2)=NTECK(ITSIDE(2,IS),ITET)
            IP(3)=NTECK(ITSIDE(3,IS),ITET)

            CUR => COORTET(IP(1))%PTET
            WHLOOP:DO WHILE (ASSOCIATED(CUR))
              JTET = CUR%NOTET
              IF (JTET /= ITET) THEN  ! OMIT TETRAHEDRON ITET
                JSLOOP:DO JS=1,4                      ! CHECK ALL SIDES
                  IF (NTBAR(JS,JTET) == 0) THEN
                    JP(1)=NTECK(ITSIDE(1,JS),JTET)
                    JP(2)=NTECK(ITSIDE(2,JS),JTET)
                    JP(3)=NTECK(ITSIDE(3,JS),JTET)
                     iloop:do i=1,3
                      jloop:do j=i,3
                        if (ip(i) == jp(j)) then
                          ip1=jp(j)
                          jp(j)=jp(i)
                          jp(i)=ip1
                          cycle iloop
                        endif
                      end do jloop
                      if (.true.) cycle jsloop
                     end do iloop
                    NTBAR(IS,ITET) = JTET ! NEIGHBOR FOUND
                    NTSEITE(IS,ITET) = JS
                    NTBAR(JS,JTET) = ITET
                    NTSEITE(JS,JTET) = IS
                    EXIT WHLOOP
                  END IF
                END DO JSLOOP ! JS
              END IF
              CUR => CUR%NEXT_TET
            END DO WHLOOP ! WHILE
          END IF
        END DO
      END DO

!      DO IC=1,NCOOR
!        CUR => COORTET(IC)%PTET
!        DO WHILE (ASSOCIATED(CUR))
!           CUR2 => CUR
!           CUR => CUR%NEXT_TET
!           DEALLOCATE (CUR2)
!        END DO
!        NULLIFY(COORTET(IC)%PTET)
!      END DO
C Allocation took already place in make_tetra.f, but output within
C make_tetra.f is difficult. Therefore this location is best for output
C of the allocated memory:
      WRITE (IUNMEM,'(A,T25,I15)') ' Nachbar-Liste ',MCLSTR*8

      RETURN
      END
