;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;(a)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;red.pro modified slightly by me to be a function usable in this
;code.  Used to calculate the luminosity distance of the sources at
;the 10 selected values of z -> Set to take the redshift, Omega_m0,
;Omega_lambda0, H_0 / 100, and a boolean verbose parameter as input
;parameters.  Set to output the luminosity distance.

function red_as_funct, z, omega0, omegalambda, h100, verbose

    common cosmology, olh

    forward_function asinh, cube, cuberoot, epeebles, dangular, dangulardiff, $
     dlosfunc, dcomovinglos, dcomovingtransverse, dhubble, thubble, $
     dluminosity,dmodulus, dvcomoving, vcomoving, $
     agefunc, getage, getredshift

; if no cosmology has been defined then define it
    
    if n_elements(olh) eq 0L then begin

       resolve_routine, 'cosmology_routines', /compile_full_file
       
       olh = {name: 'Cosmology Parameters', $ ; cosmology structure
              omega0: omega0, $
              omegalambda: omegalambda, $
              h100: h100}

    endif else begin

       IF n_elements(omega0)      NE 0 THEN olh.omega0 = omega0
       IF n_elements(omegalambda) NE 0 THEN olh.omegalambda = omegalambda
       IF n_elements(h100)        NE 0 THEN olh.h100 = h100       

    endelse 

; default cosmology
    
       
    if verbose eq 1 then begin
       print, 'Omega Matter = '+$
        strcompress(string(olh.omega0,format='(F7.4)'),/remove)
       print, 'Omega Lambda = '+$
        strcompress(string(olh.omegalambda,format='(F7.4)'),/remove)
       print, 'H_0 / 100    = '+$
        strcompress(string(olh.h100,format='(F7.4)'),/remove)
    endif 

    nz = n_elements(z)
    if nz ne 0L then for i = 0L, nz-1L do begin

       agez0 = getage(0.0)
;      print, 'Cosmology (OLH) = '
       print, 'Age (z=0.000) = '+strcompress(agez0,/remove)+' Gyr'

       if z[i] gt 0.0 then begin

          agez = getage(z[i])
          dmod = dmodulus(z[i])
          dang = dangular(z[i],/mpc)
          dlum = dluminosity(z[i],/mpc)
          xscl = dangular(z[i],/kpc)/206265.0D

          print, 'Age (z='+strcompress(string(z,format='(F5.3)'),/remove)+') = '+$
            strcompress(agez,/remove)+' Gyr'
          print, 'Lookback time = '+strcompress(agez0-agez,/remove)+' Gyr'
          print, 'DModulus      = '+strcompress(dmod,/remove)+' mag'
          print, 'DAngular      = '+strcompress(dang,/remove)+' Mpc'
          print, 'DLuminosity   = '+strcompress(dlum,/remove)+' Mpc'
          print, 'Scale         = '+strcompress(xscl,/remove)+' kpc/arcsec'

       endif
       
    endfor

return, dlum 
end  


pro homework5

;Declaring constants
c=2.998e5 ;km/s
H_0=70 ;km/s/Mpc
q_0_bench=0.5*0.3-0.7   ;(1/2)*Sum[Omega_w0*(1+3*w)] = Omega_r0+(1/2)*Omega_m0-Omega_lambda0 -> Omega_r0 neglected due to its very small value (8.4x10^-5)
                        ;w_matter=0, Omega_m0=0.3
                        ;w_lambda=-1, Omega_lambda0=0.7
                        ;w_rad=1/3, Omega_r0=8.4x10^-5 
q_0_matter=0.5*1.0 ;flat, matter-only universe -> w_matter=0 Omega_m0=1.0

;Generating 5000 random redshift values between 0 and 1
z=randomu(seed, 5000)


;Declaring my arrays of luminosity distance (5000 elements each) for
;the benchmark universe and for the flat, matter-only universe
dL_t0_bench=dblarr(n_elements(z))
dL_t0_matter=dblarr(n_elements(z))

;Using the generated redshift values to calculate the luminosity
;distances for each universe
for i=0, n_elements(dL_t0_bench)-1, 1 do dL_t0_bench[i]=(c/H_0)*z[i]*(1-((1+q_0_bench)/2)*z[i])*(1+z[i])
for j=0, n_elements(dL_t0_matter)-1, 1 do dL_t0_matter[j]=(c/H_0)*z[j]*(1-((1+q_0_matter)/2)*z[j])*(1+z[j])

;Plotting the luminosity distance as a function of redshift for both universes
cgwindow, 'cgplot', z, dL_t0_bench, symcolor=['dodger blue'], psym=1, symsize=0.5, title='Luminosity Distance as a Function of Redshift', xtitle='Redshift', ytitle='Luminosity Distance (Mpc)'

cgplot, z, dL_t0_matter, symcolor=['orange'], psym=1, symsize=0.5, /overplot, /addcmd


;Picking 10 values of z spanning the range 0<z<1 -> 0.05 to 0.95 in increments of 0.1
z_pick=dindgen(10)/10+0.05

;To calculate the luminosity distance to sources with these values of
;z, I am going to use red.pro (modified slightly here by me to make it
;a function usable in this code) -> function defined at the top of the script.

lum_dist_bench=dblarr(10)
lum_dist_matter=dblarr(10)

print, ' '

;Calculating the luminosity distances for the benchmark model ->
;Omega0=0.3, Omegalambda=0.7, h100=H_0/100=0.7
;red_as_funct(z, omega0, omegalambda, h100, verbose)
for i=0, n_elements(lum_dist_bench)-1, 1 do begin
lum_dist_bench[i]=red_as_funct(z_pick[i], 0.3, 0.7, 0.7, 1)
print, ' '
endfor

;Calcluting the luminosity distances for the matter-only model ->
;Omega0=1.0, Omegalambda=0, h100=H_0/100=0.7
for j=0, n_elements(lum_dist_matter)-1, 1 do begin
lum_dist_matter[j]=red_as_funct(z_pick[j], 1.0, 0.0, 0.7, 1)
print, ' '
endfor

;Plotting the two groups of 10 redshifts each.
cgplot, z_pick, lum_dist_bench, psym=7, symcolor='red', thick=2, /addcmd, /overplot
cgplot, z_pick, lum_dist_matter, psym=7, symcolor='charcoal', thick=2, /addcmd, /overplot

;Creating a legend to differentiate the two curves and the two groups
;of 10 redshifts.
al_legend, ['Benchmark Universe', 'Flat, Matter-Only Universe', 'Benchmark Universe', 'Flat, Matter-Only Universe'], colors=['dodger blue', 'orange', 'red', 'charcoal'], psym=[1,1,7,7], /window

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;(b)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;BENCHMARK COMPARISON

dL_t0_bench_chz=dblarr(n_elements(dL_t0_bench))
comparison_bench=dblarr(n_elements(z))
z_range_bench_new=dblarr(n_elements(comparison_bench))

for i=0, n_elements(dL_t0_bench_chz)-1, 1 do begin
dL_t0_bench_chz[i]=(c*z[i])/H_0  ;Calculating the luminosity distance given the new approximation
comparison_bench[i]=abs(dL_t0_bench[i]-dL_t0_bench_chz[i])/((dL_t0_bench[i]+dL_t0_bench_chz[i])/2)    ;Calculating the percent difference between this appraximation and the one in (a)
if comparison_bench[i] le 0.1 then z_range_bench_new[i]=z[i] else z_range_bench_new[i]=0.0
endfor

;Giving the range of z values where both approximations agree to
;within 10%
print, ' '
print, 'These two values agree to within 10% in the range: 0.0 < z < ', strcompress(string(max(z_range_bench_new)), /remove_all)
print, ' '

;Plotting the two approximations
cgwindow, 'cgplot', z, dL_t0_bench, psym=1, symsize=0.5, symcolor='dodger blue', title='Luminosity Distance as a Function of Redshift for the Two Approximations', ytitle='Luminosity Distance (Mpc)', xtitle='Redshift'

cgplot, z, dL_t0_bench_chz, psym=1, symsize=0.5, symcolor='black', /overplot, /addcmd

al_legend, ['Benchmark Taylor Expansion Approx.', 'Benchmark Simple Approx.'], psym=[1,1], color=['dodger blue', 'black'], /window


;MATTER-ONLY 

dL_t0_matter_chz=dblarr(n_elements(dL_t0_matter))
comparison_matter=dblarr(n_elements(z))
z_range_matter_new=dblarr(n_elements(comparison_matter))

for i=0, n_elements(dL_t0_matter_chz)-1, 1 do begin
dL_t0_matter_chz[i]=(c*z[i])/H_0  ;Calculating the luminosity distance given the new approximation
comparison_matter[i]=abs(dL_t0_matter[i]-dL_t0_matter_chz[i])/((dL_t0_matter[i]+dL_t0_matter_chz[i])/2)    ;Calculating the percent difference between this appraximation and the one in (a)
if comparison_matter[i] le 0.1 then z_range_matter_new[i]=z[i] else z_range_matter_new[i]=0.0
endfor

;Giving the range of z values where both approximations agree to
;within 10%
print, ' '
print, 'These two values agree to within 10% in the range: 0.0 < z < ', strcompress(string(max(z_range_matter_new)), /remove_all)
print, ' '

;Plotting the two approximations
cgwindow, 'cgplot', z, dL_t0_matter, psym=1, symsize=0.5, symcolor='red', title='Luminosity Distance as a Function of Redshift for the Two Approximations', ytitle='Luminosity Distance (Mpc)', xtitle='Redshift'

cgplot, z, dL_t0_matter_chz, psym=1, symsize=0.5, symcolor='black', /overplot, /addcmd

al_legend, ['Matter-Only Taylor Expansion Approx.', 'Matter-Only Simple Approx.'], psym=[1,1], color=['red', 'black'], /window

dL_bench=dblarr(n_elements(z))
dL_matter=dblarr(n_elements(z))

;Calculating the luminosity distance for the generated z values
;between 0 and 1 using the slightly modified red.pro
for i=0, n_elements(dL_bench)-1, 1 do dL_bench[i]=red_as_funct(z[i], 0.3, 0.7, 0.7, 0)

for i=0, n_elements(dL_matter)-1, 1 do dL_matter[i]=red_as_funct(z[i], 1.0, 0.0, 0.7, 0)

dL_bench_temp=dL_bench
dL_matter_temp=dL_matter

;Determing the array indices that correspond to luminosity distances
;greater than 100 Mpc

dL_bench_100=where(dL_bench gt 100.0)
dL_matter_100=where(dL_matter gt 100.0)

;Removing any array elements with values greater than 100 Mpc
remove, dL_bench_100, dL_bench_temp
remove, dL_matter_100, dL_matter_temp

;Determining the array index corresponding to the luminosity distance
;closest to 100 Mpc
closest_match_bench_index=where(dL_bench eq max(dL_bench_temp))
closest_match_matter_index=where(dL_matter eq max(dL_matter_temp))

print, 'The redshift at 100 Mpc in the benchmark model is: ', strcompress(string(z[closest_match_bench_index]), /remove_all)
print, ' ' 
print, 'The redshift at 100 Mpc in the matter-only model is: ', strcompress(string(z[closest_match_matter_index]), /remove_all)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;(e)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Array is 8 elements long since 5 SN 1a were excluded due to redshift
;cuts or exclusions in Hamuy et al. 1995
log_cz=dblarr(8)
b_max=dblarr(8)
v_max=dblarr(8)
m_15b=dblarr(8)
mag_max_b=dblarr(8)
mag_max_v=dblarr(8)
m_M_b=dblarr(8)
m_M_v=dblarr(8)

;Inputing the SN parameter values

log_cz[0]=4.079
log_cz[1]=4.176
log_cz[2]=4.227
log_cz[3]=3.971
log_cz[4]=4.124
log_cz[5]=3.879
log_cz[6]=4.351
log_cz[7]=4.485

b_max[0]=17.39
b_max[1]=17.88
b_max[2]=18.03
b_max[3]=16.32
b_max[4]=17.65
b_max[5]=16.11
b_max[6]=18.54
b_max[7]=19.40

v_max[0]=17.27
v_max[1]=17.83
v_max[2]=17.91
v_max[3]=16.33
v_max[4]=17.53
v_max[5]=16.13
v_max[6]=18.46
v_max[7]=19.34

m_15b[0]=1.21
m_15b[1]=1.56
m_15b[2]=1.11
m_15b[3]=1.03
m_15b[4]=1.35
m_15b[5]=1.17
m_15b[6]=1.18
m_15b[7]=1.32

;Calculating the maximum magnitudes of the SN
;From Phillips et al. 1993, M_max=a+b*m_15b

a_b=-21.726
b_b=2.698
a_v=-20.883
b_v=1.949

for i=0, n_elements(mag_max_b)-1, 1 do mag_max_b[i]=a_b+b_b*m_15b[i]
for j=0, n_elements(mag_max_v)-1, 1 do mag_max_v[j]=a_v+b_v*m_15b[j]

;Calculating m-M

for i=0, n_elements(m_M_b)-1, 1 do m_M_b[i]=b_max[i]-mag_max_b[i]
for j=0, n_elements(m_M_v)-1, 1 do m_M_v[j]=v_max[j]-mag_max_v[j]

;m-M=5*log_cz-5*log(H_0)+25
for i=0, n_elements(log_cz)-1, 1 do log_cz[i]=5*log_cz[i]

;Plotting m-M vs. log(cz)

cgwindow, 'cgplot', log_cz, m_M_b, psym=1, symcolor=['green'], thick=2, title='Hubble Diagram in the B-Band', ytitle='m-M', xtitle='log_10(cz)', yrange=[34, 38]

;Calculating the slope and y-intercept of the best fit line
best_fit_b=dblarr(2)

best_fit_b=linfit(log_cz, m_M_b)

;linfit returns [y-int, slope]
;Calculating H_0 as described in part d
H_0_b=10^((best_fit_b[0]-25)/(-5))

print, ' '
print, 'H_0 for the B-band is: ', strcompress(string(H_0_b), /remove_all)

;Creating the x values for the best fit equation from which the y
;values are calculated.
x_H_0=dindgen(1000)/1000*4+19
y_H_0_b=best_fit_b[1]*x_H_0+best_fit_b[0]

;Overplotting the best fit line
cgplot, x_H_0, y_H_0_b, psym=1, symcolor=['red'], symsize=0.5,  /overplot, /addcmd

;Calculating the slope and y-intercept of the best fit line
best_fit_v=dblarr(2)

best_fit_v=linfit(log_cz, m_M_v)

;linfit returns [y-int, slope]
;Calculating H_0 as described in part d
H_0_v=10^((best_fit_v[0]-25)/(-5))

print, ' '
print, 'H_0 for the V-band is: ', strcompress(string(H_0_v), /remove_all)

print, ' '

;Plotting m-M vs. log(cz)
cgwindow, 'cgplot', log_cz, m_M_v, psym=1, symcolor=['brown'], thick=2, title='Hubble Diagram in the V-Band', ytitle='m-M', xtitle='log_10(cz)', yrange=[34,38]


;Calculating the y values with the best fit equation and the
;previously created x values
y_H_0_v=best_fit_v[1]*x_H_0+best_fit_v[0]

;Overplotting the best fit line

cgplot, x_H_0, y_H_0_v, psym=1, symcolor=['red'], symsize=0.5, /overplot, /addcmd

end

