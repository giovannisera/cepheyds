!MODULI
MODULE strumenti
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE w_mean(x, sigma, n, media, media_err)
        IMPLICIT NONE
    
        INTEGER :: n
        REAL*8 :: x(n), sigma(n)
    
        REAL*8 :: num, denom
        REAL*8 :: w
        REAL*8, INTENT(OUT) :: media, media_err
    
        INTEGER :: k
    
        num = 0.d0
        denom = 0.d0
        DO k=1, n
            w = 1.d0 /(sigma(k))**2.d0 !peso
            num = num + x(k) *w 
            denom = denom + w 
        ENDDO
        media = num /denom
        media_err = 1.d0 /SQRT(denom)
    END SUBROUTINE w_mean

    SUBROUTINE sorting(x, ndati, y, z)
        IMPLICIT NONE
    
        INTEGER :: ndati, posizione
        INTEGER :: i, j
    
        REAL*8 :: minimo, yp, zp
        REAL*8 :: x(ndati)
        REAL*8, OPTIONAL :: y(ndati), z(ndati)
        
        DO i=1, ndati-1
            minimo = x(i)
            posizione = i
            DO j=i+1, ndati
                IF (x(j) < minimo) THEN
                    minimo = x(j)
                    posizione = j
                ENDIF
            ENDDO
            x(posizione) = x(i)
            x(i) = minimo
            IF (PRESENT(y)) THEN
                yp = y(posizione)
                y(posizione) = y(i)
                y(i) = yp
                IF (PRESENT(z)) THEN
                    zp = z(posizione)
                    z(posizione) = z(i)
                    z(i) = zp
                ENDIF
            ENDIF
        ENDDO
    END SUBROUTINE sorting
END MODULE strumenti

MODULE fase_1
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE jordan(a, c, neq, x)
        IMPLICIT NONE

        INTEGER :: neq
        INTEGER :: i, j, k 

        REAL*8 :: fakt, aux
        REAL*8 :: a(neq,neq), c(neq)
        REAL*8, INTENT(OUT) :: x(neq)
    
        DO i=1, neq 
            aux = a(i,i)
            DO j=1, neq 
                a(i,j) = a(i,j) / aux
            ENDDO
            c(i) = c(i) / aux
            DO j=1, neq 
                IF (i /= j) THEN
                    fakt = a(j,i) / a(i,i)
                    DO k=1, neq 
                        a(j,k) = a(j,k) - a(i,k)*fakt
                    ENDDO
                    c(j) = c(j) - fakt*c(i)
                ENDIF
            ENDDO
        ENDDO

        x = c
    END SUBROUTINE jordan

    SUBROUTINE fit_lin(x, y, ndati, nv, coeff)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ndati, nv
        INTEGER :: i, j, k

        REAL*8 :: summa
        REAL*8 :: a(nv,nv), c(nv) 
        REAL*8, INTENT(IN) :: x(ndati), y(ndati)
        REAL*8, INTENT(OUT) :: coeff(nv)

        DO i=1, nv
            DO j=1, nv
                SUMMA = 0.d0
                DO k=1, ndati
                    summa = summa + x(k)**(i+j-2)
                ENDDO
                a(i,j) = summa
            ENDDO
        ENDDO
    
        DO i=1, nv
            summa = 0.d0
            DO k=1, ndati
                summa = summa + y(k)*(x(k)**(i-1))
            ENDDO
            c(i) = summa
        ENDDO

        CALL jordan(a, c, nv, coeff)
    END SUBROUTINE fit_lin
END MODULE fase_1

MODULE fase_2
    USE fase_1
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE inter(x1, x2, dx2, y, n, out)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        INTEGER :: i

        REAL*8, INTENT(IN) :: x1(n), x2(n), dx2(n), y
        REAL*8, INTENT(OUT) :: out

        DO i=2, n-1
            IF (y >= x1(i-1) .and. y <= x1(i)) EXIT
        ENDDO

        out = (dx2(i-1) / (6.d0* (x1(i) - x1(i-1)))) * (x1(i) - y)**3 +&
        (dx2(i) / (6.d0* (x1(i) - x1(i-1)))) * (y - x1(i-1))**3 +&
            ((x2(i-1) / (x1(i) - x1(i-1))) - dx2(i-1) * (x1(i)-x1(i-1)) /6.d0) * (x1(i) - y) +&
                ((x2(i) / (x1(i) - x1(i-1))) - (dx2(i) * (x1(i) - x1(i-1))) /6.d0) * (y - x1(i-1)) 
    END SUBROUTINE inter
 
    SUBROUTINE spline_cubica(x1, x2, n, dx2)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        INTEGER :: i, k

        REAL*8 :: c(n), m(n,n)
        REAL*8, INTENT(IN) :: x1(n), x2(n)
        REAL*8, INTENT(OUT) :: dx2(n)
    
        m = 0.d0
	    m(1,1) = 1.d0
	    m(n,n) = 1.d0
	    DO i=2, n-1
		    m(i,i-1) = x1(i) - x1(i-1) !e
		    m(i,i) = 2.d0* (x1(i+1) - x1(i-1)) !r
		    m(i,i+1) = x1(i+1) - x1(i) !g
        ENDDO

        c(1) = 0.d0
	    c(n) = 0.d0
	    DO i=2, n-1
		    c(i) = (6.d0/ (x1(i+1) - x1(i))) * (x2(i+1) - x2(i)) +&
                    (6.d0/ (x1(i) - x1(i-1))) * (x2(i-1) - x2(i))
        ENDDO

        CALL jordan(m, c, n, dx2)
    END SUBROUTINE spline_cubica

    SUBROUTINE chiquadro(periodi_test, chi, time, appmag, err, len, n, dx2)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: n, len
        INTEGER :: i, j, k, count

        REAL*8 :: summa, out, steps(2)
        REAL*8, INTENT(IN) :: time(len), appmag(len), err(len), dx2(len)
        REAL*8, INTENT(OUT) :: periodi_test(n), chi(n)  
        
        DO i=1, n !scelgo un periodo da testare
            summa = 0.d0
            count = 0 
            !associo ad ogni periodo di test un chi-quadro iterando sull'intervallo temporale
            DO j=1, len 
                steps(1) = time(j) + periodi_test(i)
                steps(2) = time(j) - periodi_test(i)
                DO k=1, 2
                    DO WHILE ((steps(k) >= time(1)) .and. (steps(k) <= time(len))) 
                        CALL inter (time, appmag, dx2, steps(k), len, out)
                        summa = summa + ((out - appmag(j)) / err(j))**2
                        SELECT CASE(k)
                        CASE(1)
                            steps(k) = steps(k) + periodi_test(i)
                        CASE(2)
                            steps(k) = steps(k) - periodi_test(i)
                        END SELECT
                        count = count + 1
                    ENDDO
                ENDDO
            ENDDO
            chi(i) = summa / (count - 1) !chiquadro ridotto
        ENDDO    
    END SUBROUTINE chiquadro
END MODULE fase_2

MODULE fase_cosmologia
    IMPLICIT NONE
    CONTAINS
        SUBROUTINE quadratura (integranda, a, b, OM, OL, risultato)
            IMPLICIT NONE
            
            INTEGER, PARAMETER :: npunti = 4
            INTEGER :: i

            REAL*8, INTENT(IN) :: a, b
            REAL*8 :: OM, OL, summa
            REAL*8, INTENT(OUT) :: risultato
            REAL*8, ALLOCATABLE :: c(:), xd(:), x(:)
            REAL*8, EXTERNAL :: integranda

            ALLOCATE(c(0:npunti-1), xd(0:npunti-1), x(0:npunti-1))

            c(0) = (18.d0 + sqrt(30.d0)) /36.d0
            c(1) = c(0)
            c(2) = (18.d0 - sqrt(30.d0)) /36.d0
            c(3) = c(2)
            xd(0) = SQRT((3.d0 /7.d0) - (2.d0 /7.d0) * SQRT(6.d0 /5.d0))
            xd(1) = -xd(0)
            xd(2) = SQRT((3.d0 /7.d0) + (2.d0 /7.d0) *SQRT(6.d0 /5.d0))
            xd(3) = -xd(2)

            DO i=0, npunti-1
                x(i) = ((b+a) + (b-a) *xd(i)) /2.d0 !cambio di variabile per calcolare l'integrale tra -1 e +1
            ENDDO 
     
            summa = 0.d0
     
            DO i=0, npunti-1
                summa = summa + c(i) *integranda(OM, OL, x(i))
            ENDDO
            risultato = summa *(b-a) /2.d0 ! moltiplico per il contributo dovuto al differenziale del cambio di variabile
        END SUBROUTINE quadratura
END MODULE fase_cosmologia

!PROGAMMA
PROGRAM main    
    USE strumenti
    USE fase_1
    USE fase_2
    USE fase_cosmologia
    IMPLICIT NONE

    CHARACTER(LEN=1) :: letters(3)
    CHARACTER(LEN=4) :: galaxies(10)
    CHARACTER(LEN=20), EXTERNAL :: str_int, str_real
    CHARACTER(LEN=7) :: nulla

    INTEGER :: len, start, n, nn, mm
    INTEGER, PARAMETER :: nv = 2, npoints = 1000
    INTEGER :: i, j, k, h, g, scorri, l, count, count2 !per iterazioni
    
    REAL*8 :: p_min, p_max, out, out_sx, out_dx, psx_min, psx_max, pdx_min, pdx_max,&
                summa, summa2, pe, w, hubble_constant,&
                    err_hubble_constant, OM, OL, t_0, t_1, t_2, t, alpha

    REAL*8 :: coeff(nv), x(0:npoints-1), y(0:npoints-1), periodo(30), err_periodi(30), magnitudini_apparenti(30),&
                errori_magnitudini_apparenti(30), magnitudini_assolute(30), errori_magnitudini_assolute(30),&
                    modulo_di_distanza(30), errori_modulo_dist(30), distanze(30), err_distanze(30), dist(10),& 
                        err_dist(10), vel_rec(10), err_vel_rec(10), H0(10), err_H0(10), dist_gal(3), err_dist_gal(3)

    REAL*8, ALLOCATABLE :: P(:), M(:), time(:), appmag(:), err_appmag(:),&
                            dx2(:), periodi_test(:), chi(:), lower(:), p_sx(:),&
                            chi_sx(:), p_dx(:), chi_dx(:), media(:), RMS(:), spline(:)

    REAL*8, EXTERNAL :: integranda1, integranda2

    PRINT*, "#Stima della costante di Hubble, partendo dalle curve di luce" 
    PRINT*, "di 30 stelle cefeidi osservate in 10 diverse galassie#"
    PRINT*, ""
    PRINT*,"->INIZIO"
    PRINT*, ""

    PRINT*, "###FASE 1 in corso..."
    OPEN(1, FILE = "ceph_catalog.txt")
    len = 0
    DO 
        READ(1, *, END=100)
        len = len + 1
    ENDDO
    100 CONTINUE
    
    ALLOCATE(P(len), M(len))
    REWIND(1)
    
    !REGRESSIONE LINEARE DELLA RELAZIONE MV vs. LOG10(P)
    OPEN(2, FILE="tab_ceph.txt")
    DO i=1, len
        READ(1, *) P(i), M(i)
    ENDDO
    CLOSE(1)
    CALL sorting(P, len, M)
    DO i=1, len
        WRITE(2, *) P(i), M(i)
    ENDDO
    CLOSE(2)

    CALL fit_lin(LOG10(P), M, len, nv, coeff)
    PRINT*, "I coefficienti della relazione M_V vs. logP sono stati trovati"
    PRINT*, "coeff"//trim(str_int(1))//" = "//trim(str_real(coeff(1)))//&
                    ", coeff"//trim(str_int(2))//" = "//trim(str_real(coeff(2)))
    PRINT*, "###FINE FASE 1"
    PRINT*, ""

    PRINT*, "###FASE 2 in corso..."

    OPEN(7, FILE = "gal_vel.txt")
    !"scarico" sulla variabile "nulla" i nomi delle galassie, non necessari
    DO i=1, 10
        READ(7, *) nulla, vel_rec(i), err_vel_rec(i) 
    ENDDO
    CLOSE(7)
    galaxies = (/ "0925", "1365", "1425", "2403", "2541", "3198", "3621",&
                    "4535", "4536", "4639"/)
    letters = (/"a", "b", "c"/)
    scorri = 1
    DO i=1, 10
        DO j=1, 3
            !LETTURA ED ORDINAMENTO IN BASE AL TEMPO DEL FILE DI OGNI CEFEIDE
            OPEN(1, FILE = "ceph_NGC"//galaxies(i)&
                            //letters(j)//".txt")
            len = 0
            DO 
                READ(1, *, END=101)
                len = len + 1
            ENDDO
            101 CONTINUE
            ALLOCATE(time(len), appmag(len), err_appmag(len), dx2(len))
            REWIND(1)
            OPEN(2, FILE="tab_ceph"//galaxies(i)&
                                //letters(j)//".txt")
            DO k=1, len
                READ(1, *) time(k), appmag(k), err_appmag(k)
            ENDDO
            CLOSE(1)
            CALL sorting(time, len, appmag, err_appmag)
            DO k=1, len
                WRITE(2, *) time(k), appmag(k), err_appmag(k)
            ENDDO
            CLOSE(2)

            !SPLINE CUBICA DELLA CURVA DI LUCE
            CALL spline_cubica(time, appmag, len, dx2)
            OPEN(3, FILE="spline_ceph"//galaxies(i)&
                                //letters(j)//".txt")
            DO k=0, npoints-1
                x(k) = time(1) + ((time(len) - time(1)) / (npoints - 1)) * k
            ENDDO               
            y(0) = appmag(1)
            y(npoints-1) = appmag(len)
            DO k=1, npoints-2
                CALL inter (time, appmag, dx2, x(k), len, out)
                y(k) = out
            ENDDO
            DO k=0, npoints-1
                WRITE(3, *) x(k), y(k)
            ENDDO
            CLOSE(3)

            !CALCOLO DEL PERIODO CON ERRORE
            !cerco la minima distanza temporale lower(1) nell'array che contiene le date
            start = 0
            !inserisco in un array tutte le possibili combinazioni delle differenze tra i tempi
            !lunghezza dell'array (combinazioni di len elementi senza ripetione in classe due)
            DO k=1, len
                start = start + (len-k)
            ENDDO
            ALLOCATE(lower(start))
            k = 1
            DO h=1, len
                DO g=h+1, len
                    lower(k) = ABS(time(h) - time(g))
                    k = k + 1
                ENDDO
            ENDDO
            CALL sorting(lower, start)
            p_min = 3.d0 * lower(1) !periodo minimo da testare
            p_max = (time(len) - time(1)) / 2.d0 !periodo massimo da testare
            !diviso l'asse temporale in n intervallini lunghi 0.02 giorni
            n = (p_max - p_min) / 0.02d0
            ALLOCATE(periodi_test(0:n-1), chi(0:n-1))
            !creo la griglia temporale
            DO k=0, n-1
                periodi_test(k) = p_min + 0.02d0 * k
            ENDDO
            OPEN(4, FILE="chiquadro_ceph"//galaxies(i)&
                                //letters(j)//".txt")
            CALL chiquadro(periodi_test, chi, time, appmag, err_appmag, len, n, dx2)
            DO k=0, n-1
                WRITE(4, *) periodi_test(k), chi(k)
            ENDDO
            CLOSE(4)
            CALL sorting(chi, n, periodi_test)
            !errore sul periodo: spline inversa nella curva del chi-quadro, interpolo a sx e a dx del periodo minimo
            !SX (stimo l'errore essere minore di 0.2 giorni dopo aver osservato i grafici del chi-quadro)
	        psx_min = periodi_test(0) - 0.2d0 
	        psx_max = periodi_test(0)
	        nn = (psx_max - psx_min) / 0.02d0
	        ALLOCATE(p_sx(0:nn-1), chi_sx(0:nn-1))
	        DO k=0, nn-1
		        p_sx(k) = psx_min + 0.02d0 * k
	        ENDDO
	        CALL chiquadro(p_sx, chi_sx, time, appmag, err_appmag, len, nn, dx2)
            CALL sorting(chi_sx, size(chi_sx), p_sx)
            CALL spline_cubica(chi_sx, p_sx, size(chi_sx), dx2)
            !1 sigma a sinistra
	        CALL inter(chi_sx, p_sx, dx2, chi(0) + 1.d0, size(chi_sx), out_sx)
	        !DX 
	        pdx_min = periodi_test(0) 
	        pdx_max = periodi_test(0) + 0.2d0 
	        nn = (pdx_max - pdx_min) / 0.02d0
	        ALLOCATE(p_dx(0:nn-1), chi_dx(0:nn-1))
	        DO k=0, nn-1
		        p_dx(k) = pdx_min + 0.02d0 * k
	        ENDDO
            !devo ricavarmi di nuovo dx2 della curva di luce per calcolare il chi-quadro
            CALL spline_cubica(time, appmag, len, dx2) 
	        CALL chiquadro(p_dx, chi_dx, time, appmag, err_appmag, len, nn, dx2)
            CALL sorting(chi_dx, size(chi_dx), p_dx)
            CALL spline_cubica(chi_dx, p_dx, size(chi_dx), dx2)
            !1 sigma a destra
	        CALL inter(chi_dx, p_dx, dx2, chi(0) + 1.d0, size(chi_dx), out_dx)
	        !l'errore sul periodo è il massimo tra i due errori trovati
	        IF (out_dx >= out_sx) THEN
		        err_periodi(scorri) = ABS(out_dx - periodi_test(0))
	        ELSE
		        err_periodi(scorri) = ABS(out_sx - periodi_test(0))
	        ENDIF
            periodo(scorri) = periodi_test(0)
 
            !MAGNITUDINE APPARENTE MEDIA CON ERRORE ASSOCIATO
            CALL spline_cubica(time, appmag, len, dx2) 
            pe = periodi_test(0) !periodo
            n = INT((time(len) - time(1)) /(pe /100.d0)) + 1 !numero di dati da interpolare
            mm = INT((time(len) - time(1)) /pe) + 1 !numero di periodi coperti dai dati
            ALLOCATE(spline(n))
            ALLOCATE(media(mm), RMS(mm))
            t = time(1) !istante di tempo iniziale
            count = 0 !tengo conto dei dati presenti in ciascun periodo
            count2 = 1 !tengo conto delle interpolazioni su ciascun periodo
            h = 1 !per scorrere su media e RMS
            k = 1 !per scorrere su spline
            g = 1 !tengo conto dei periodi coperti
            l = 1 !per scorrere su spline #2
            summa = 0.d0 !per il calcolo della media
            summa2 = 0.d0 !per il calolo dell'errore
            DO WHILE (t <= time(len))
                CALL inter(time, appmag, dx2, t, len, out)
                spline(k) = out
                summa = summa + spline(k)
                IF ((t - time(1)) >= h *pe .or. t + (pe /100.d0) > time(len)) THEN
                    media(h) = summa /count
                    DO l=k, count2, -1
                        summa2 = summa2 + (spline(l) - media(h))**2
                    ENDDO
                    RMS(h) = SQRT(summa2 /count)
                    summa = 0.d0
                    summa2 = 0.d0
                    count = 0
                    count2 = k+1
                    h = h+1
                    g = g+1
                ENDIF
                k = k+1
                count = count+1
                t = t + pe/100.d0
            ENDDO
            CALL w_mean(media, RMS, mm, magnitudini_apparenti(scorri), errori_magnitudini_apparenti(scorri))

            !MAGNITUDINE ASSOLUTA, MODULO DI DISTANZA, DISTANZA CON ERRORE ASSOCIATO
            magnitudini_assolute(scorri) = coeff(1) + coeff(2) *LOG10(periodo(scorri))
            errori_magnitudini_assolute(scorri) = ABS(coeff(2) /(periodo(scorri) *LOG(10.d0))) *err_periodi(scorri)
            modulo_di_distanza(scorri) = magnitudini_apparenti(scorri) - magnitudini_assolute(scorri)
            errori_modulo_dist(scorri) = SQRT( errori_magnitudini_apparenti(scorri)**2 +&
                                         errori_magnitudini_assolute(scorri)**2 ) !somma in quadratura
            distanze(scorri) = (10 **(0.2d0 *modulo_di_distanza(scorri) + 1.d0))  /10**6 ![MegaParsec]
            err_distanze(scorri) = 0.2d0 *LOG(10.d0) *distanze(scorri) *errori_modulo_dist(scorri) ![MegaParsec]
            dist_gal(j) = distanze(scorri)
            err_dist_gal(j) = err_distanze(scorri)
            scorri = scorri+1

            DEALLOCATE(time, appmag, err_appmag, dx2, lower,&
                           periodi_test, chi, p_sx, chi_sx, p_dx, chi_dx,&
                                media, RMS, spline)
        ENDDO
        
        !DISTANZA E COSTANTE DI HUBBLE DELLE 10 GALASSIE
        CALL w_mean(dist_gal, err_dist_gal, 3, dist(i), err_dist(i))
        H0(i) = vel_rec(i) / dist(i)
        err_H0(i) = H0(i) *sqrt((err_vel_rec(i) /vel_rec(i))**2 + (err_dist(i) /dist(i))**2)
    ENDDO

    !SALVATAGGIO CARATTERISTICHE DELLE 30 CEFEIDI
    OPEN(5, FILE="periodi_ceph.txt")
    k = 1
    DO i=1, 10
        DO j=1, 3
            WRITE(5, *) periodo(k), err_periodi(k), magnitudini_apparenti(k),&
                            errori_magnitudini_apparenti(k), magnitudini_assolute(k),&
                                errori_magnitudini_assolute(k), modulo_di_distanza(k),&
                                    errori_modulo_dist(k), distanze(k), err_distanze(k),&
                                        "NGC"//galaxies(i)//letters(j)
            k = k + 1
        ENDDO
    ENDDO
    CLOSE(5)
    PRINT*, "Tutte le cefeidi sono state analizzate; consultare i file output e la relazione in allegato"
    PRINT*, "###FINE FASE 2"
    PRINT*, ""

    !SALVATAGGIO CARATTERISTICHE DELLE 10 GALASSIE
    PRINT*, "###FASI 3-4 in corso..."
    OPEN(100, FILE = "hubble.txt")
    DO i=1, 10
        WRITE(100, *) dist(i), err_dist(i), vel_rec(i), err_vel_rec(i),&
                     H0(i), err_H0(i), "NGC"//galaxies(i)
    ENDDO
    CLOSE(100)

    !STIMA DELLA COSTANTE DI HUBBLE CON ERRORE
    CALL w_mean(H0, err_H0, 10, hubble_constant, err_hubble_constant)
    PRINT*, "Tutte le galassie sono state analizzate; consultare i file output e la relazione in allegato"
    PRINT*, "H0="//trim(str_real(hubble_constant))
    PRINT*, "Errore associato ad H0="//trim(str_real(err_hubble_constant))
    PRINT*, "###FINE FASI 3-4"
    PRINT*, ""
    PRINT*, "###FASE COSMOLOGIA in corso..."

    !ETA' DELL'UNIVERSO E PARAMETRI DI DENSITA'
    !esprimo la costante di Hubble in [Mldyy^{-1}]
    alpha = 31536.d0 /30857.d3
    hubble_constant = hubble_constant *alpha
    OPEN(300, FILE = "1sigma.txt")
    OPEN(400, FILE = "2sigma.txt")
    OPEN(500, FILE = "3sigma.txt")
    OM = 0.d0
    DO WHILE ( OM >= 0.d0 .and. OM <= 1 )
        OL = 0.d0
        DO WHILE ( OL >= 0.d0 .and. OL <= 1 )
            CALL quadratura( integranda1, 0.d0, 1.d0, OM, OL, t_1) ! primo addendo dell'integrale
            CALL quadratura( integranda2, 0.d0, 1.d0, OM, OL, t_2) ! secondo addendo dell'integrale
            t_0 = (1.d0 /hubble_constant) * (t_1 + t_2) 
            IF (t_0 >= 13.82d0 - 0.14d0 .and. t_0 <= 13.82d0 + 0.14d0) THEN
                !1-sigma
                WRITE(300, *) OM, OL
            ELSE IF (t_0 >= 13.82d0 - 0.28d0 .and. t_0 <= 13.82d0 + 0.28d0) THEN
                !2-sigma
                WRITE(400, *) OM, OL
            ELSE IF (t_0 >= 13.82d0 - 0.42d0 .and. t_0 <= 13.82d0 + 0.42d0) THEN
                !3-sigma 
                WRITE(500, *) OM, OL
            ENDIF
            OL = OL + 0.001d0
        ENDDO
        OM = OM + 0.001d0
    ENDDO 
    CLOSE(300)
    CLOSE(400)  
    CLOSE(500)
    PRINT*, "Sono stati ricavati i modelli compatibili con Planck 2018 e H_0 di questo progetto entro 3-sigma;"
    PRINT*, "consultare i file output"
    PRINT*, "###FINE FASE COSMOLOGIA"
    PRINT*,""
    PRINT*, "AVVISO-controllare la cartella di memoria contente il programma per visualizzare gli array di dati prodotti"
    PRINT*, "AVVISO-aprire l'allegato JUPYTER per avere i plot relativi a NGC4536 sul vostro computer"
    PRINT*,""
    PRINT*,"FINE DEL PROGRAMMA<-giovanni.sera@studio.unibo.it matricola: 1030375 in data 28/08/2023"
END PROGRAM main

!DUE FUNZIONI USATE PER FORMATTARE MEGLIO LA STAMPA IN OUTPUT SUL TERMINALE
CHARACTER(LEN=20) FUNCTION str_int(k)
!   "Converte un numero intero in una stringa."
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: k

    WRITE (str_int, *) k
    str_int = ADJUSTL(str_int)

END FUNCTION str_int
CHARACTER(LEN=20) FUNCTION str_real(k)
!   "Converte un numero reale in una stringa."
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: k

    WRITE (str_real, '(F7.2)') k
    str_real = ADJUSTL(str_real)

END FUNCTION str_real

!NECESSARIE PER IL CALCOLO DELL'INTEGRALE t_0
REAL*8 FUNCTION integranda1(OM, OL, z)
! "Funzione Integranda del Primo Addendo."  
    IMPLICIT NONE

    REAL*8 ::OM, OL, z, E

    E = OM * (1.d0 + z)**3 + (1.d0 - OM - OL) * (1.d0 + z)**2 + OL
    integranda1 =  1.d0 / ( (1 + z) * SQRT(E) ) 

END FUNCTION integranda1
REAL*8 FUNCTION integranda2(OM, OL, t)
! "Funzione Integranda del Secondo Addendo."  
    IMPLICIT NONE

    REAL*8 ::OM, OL, t, E

    E = OM * (1.d0 + (1.d0 / t))**3 + (1.d0 - OM - OL) * (1.d0 + (1.d0 / t))**2 + OL ! cambio di variabile in integranda 1 z=>1/t
    integranda2 = 1.d0 / (t**2 * (1.d0 + (1.d0 / t)) * SQRT(E) ) ! funzione integranda finale che comprende anche il fattore 1/t**2
                                                                    ! derivante dal cambio di variabile

END FUNCTION integranda2

    
