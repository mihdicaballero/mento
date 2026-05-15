i nned ot work now in the shear wall summary, similar to beams.

here the input is from an excel too and is like this:

										Transverse rebar		Vertical rebar
Level	Label	Comb.	t	lw	hw	cc	Nx	Vz	My	dbh	sh	dbv	sv
			cm	m	m	mm	kN	kN	kNm	mm	cm	mm	cm
Level 1	M1	ELU 1	20	3.00	3.00	25	0	264	-172	8	20	12	15
Level 1	M1	ELU 2	20	3.00	3.00	25	0	138	-90	8	20	12	15
Level 1	M1	ELU 3	20	3.00	3.00	25	0	123	-81	8	20	12	15
Level 1	M1	ELU 4	20	3.00	3.00	25	-301	152	-234	8	20	12	15
Level 2	M1	ELU 1	20	3.00	3.00	25	-150	32.3	143	8	20	12	15
Level 2	M1	ELU 2	20	3.00	3.00	25	55.5	163	-278	8	20	12	15
Level 2	M1	ELU 3	20	3.00	3.00	25	282	19	159	8	20	12	15
Level 2	M1	ELU 4	20	3.00	3.00	25	-4.5	88.15	-97	8	20	12	15
Level 1	M2	ELU 1	20	2.00	3.00	25	-240	61.2	-38	8	20	12	15
Level 1	M2	ELU 2	20	2.00	3.00	25	-163	29	60	8	20	12	15
Level 1	M2	ELU 3	20	2.00	3.00	25	-17	47	-39	8	20	12	15
Level 1	M2	ELU 4	20	2.00	3.00	25	332	21	46.13	8	20	12	15
Level 2	M2	ELU 1	20	2.00	3.00	25	-150	32.3	143	8	20	12	15
Level 2	M2	ELU 2	20	2.00	3.00	25	55.5	163	-278	8	20	12	15
Level 2	M2	ELU 3	20	2.00	3.00	25	-163	29	60	8	20	12	15
Level 2	M2	ELU 4	20	2.00	3.00	25	55.5	163	-278	8	20	12	15


we would need to add "level" to shear wall atributes i believe. then we have same logic of input handler

the thing is that summary was made only for beams, maybe we need to refactor a beams_summary and create a shear_wall_summary (i will have punching summary too later)

most of methods form beam summary apply here, we will also do check i fthere is rebar indicated or design if its empty.

theres is a difference here: in shear wall a level can have a lsit of forces for a wall, so that wall will have a lsit of forces to be cheked or design, whereas in beams the summary was for a list of independent beam sections.

also in shear wall i'm not intereseted in the capacityu=true scensario, i just design or check with rebar.

also this is only for shear, there is no flexure in shear wall implement and eurocode is not implemented either in any case yet

i would also have to create tests with a sample dataframe for this. for this matter, create a file in _notebooks as .py so you can work in there the methods and check results.
