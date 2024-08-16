Calculations
Let the force be represented as:
$$F=ma$$
Where:
- $m$ is the mass in $\text{kg}$
- $a$ is the acceleration in $\text{m/s}^2$


The result is: $F$=5.0 kN



The result is: $I_x$=16457503 mm⁴



The result is: $I_y$=2986995 mm⁴





1.124 kip




\[
\begin{aligned}
\mathrm{mass} &= 100.00\ \mathrm{kg} \; 
\\[10pt]
a &= 9.81 \cdot \frac{ m }{ \left( s \right) ^{ 2 } }  = 9.81 \cdot \frac{ m }{ \left( s \right) ^{ 2 } } &= 9.81\ \mathrm{m} \cdot \mathrm{s}^{-2}  
\\[10pt]
\mathrm{force} &= \mathrm{mass} \cdot a  = 100.00\ \mathrm{kg} \cdot 9.81\ \mathrm{m} \cdot \mathrm{s}^{-2} &= 981.00\ \mathrm{N}  
\end{aligned}
\]



\[
\begin{aligned}
\mathrm{mass} &= 100.00\ \mathrm{kg} \; 
\\[10pt]
a &= 9.81\ \mathrm{m} \cdot \mathrm{s}^{-2} \; 
\\[10pt]
\mathrm{force} &= 981.00\ \mathrm{N} \; 
\end{aligned}
\]



\[
\begin{aligned}
\mathrm{mass} &= 100.00\ \mathrm{kg} \; 
 &a &= 9.81\ \mathrm{m} \cdot \mathrm{s}^{-2} \; 
 &\mathrm{force} &= 981.00\ \mathrm{N} \; 
\\[10pt]
\end{aligned}
\]



\[
\begin{aligned}
\mathrm{mass} &= 100 \cdot \mathrm{kg} \; \;\textrm{(Comment in line)}
\\[10pt]
a &= 9.81 \cdot \frac{ m }{ \left( s \right) ^{ 2 } } \; 
\\[10pt]
\mathrm{force} &= \mathrm{mass} \cdot a \; 
\end{aligned}
\]



\[
\begin{aligned}
\mathrm{force} &= 981.00\ \mathrm{N} \;
\end{aligned}
\]



\[
\begin{aligned}
t &= 3.75\ \mathrm{mm} \; \;\textrm{(Espesor)}
 &H &= 200.00\ \mathrm{mm} \; \;\textrm{(Altura)}
 &B &= 100.00\ \mathrm{mm} \; \;\textrm{(Ancho)}
\\[10pt]
\end{aligned}
\]



\[
\begin{aligned}
I_{x} &= B \cdot \frac{ \left( H \right) ^{ 3 } }{ 12 } - \left( B - t \right) \cdot \frac{ \left( H - t \right) ^{ 3 } }{ 12 } \\&= 100.000\ \mathrm{mm} \cdot \frac{ \left( 200.000\ \mathrm{mm} \right) ^{ 3 } }{ 12 } - \left( 100.000\ \mathrm{mm} - 3.750\ \mathrm{mm} \right) \cdot \frac{ \left( 200.000\ \mathrm{mm} - 3.750\ \mathrm{mm} \right) ^{ 3 } }{ 12 } \\&= 6042122.192\ \mathrm{mm}^{4}  \\[10pt]
\\[10pt]
I_{y} &= H \cdot \frac{ \left( B \right) ^{ 3 } }{ 12 } - \left( H - t \right) \cdot \frac{ \left( B - t \right) ^{ 3 } }{ 12 } \\&= 200.000\ \mathrm{mm} \cdot \frac{ \left( 100.000\ \mathrm{mm} \right) ^{ 3 } }{ 12 } - \left( 200.000\ \mathrm{mm} - 3.750\ \mathrm{mm} \right) \cdot \frac{ \left( 100.000\ \mathrm{mm} - 3.750\ \mathrm{mm} \right) ^{ 3 } }{ 12 } \\&= 2084212.036\ \mathrm{mm}^{4}  \\[10pt]
\end{aligned}
\]



    ---------------------------------------------------------------------------

    TypeError                                 Traceback (most recent call last)

    Cell In[44], line 42
         33     section = Beam(
         34         name="V-40x50",
         35         concrete=concrete,
       (...)
         38         depth=500 * mm,  # type: ignore
         39     )
         41 if __name__ == "__main__":
    ---> 42     main()
    

    Cell In[44], line 32, in main()
         29 def main():
         30     # Ejemplo de uso
         31     concrete=create_concrete(name="H30",f_c=30*MPa, design_code="ACI 318-19") # type: ignore
    ---> 32     steelBar=SteelBar(name="ADN 420", f_y=420*MPa) # type: ignore
         33     section = Beam(
         34         name="V-40x50",
         35         concrete=concrete,
       (...)
         38         depth=500 * mm,  # type: ignore
         39     )
    

    TypeError: SteelBar.__init__() got an unexpected keyword argument 'f_y'

