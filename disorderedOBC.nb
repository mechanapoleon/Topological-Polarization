(* PARAMETERS *)
charge     = 1.0;  (* charge *)
aa         = 10;   (* lattice parameter *)
\[Epsilon] = 1;    (* energy per site *)
t          = -1.0; (* hopping energy *)


(* FUNCTIONS *)

(* Hamiltonian with open boundary conditions *)
hamchainOBCdisorder[sites_] := SparseArray[Which[sites == 1, {{i_, i_} -> +RandomReal[{0.5, 1}]}, sites == 2, {{1, 1} -> +RandomReal[{0.5, 1}], {2, 2} -> -RandomReal[{0.5, 1}], {1, sites} -> t, {sites, 1} -> t}, sites > 2, {{i_?OddQ, i_?OddQ} -> RandomReal[{0.5, 1}], {i_?EvenQ, i_?EvenQ} -> -RandomReal[{0.5, 1}], {i_, j_} /; Or[i == (j + 1), i == (j - 1), j == (i + 1), j == (i - 1)] -> t}], {sites, sites}];

(* Energy Eigenvalues *)
energy[hamiltonian_] := Eigensystem[hamiltonian][[1]];

(* Energy Eigenvectors *)
vects[hamiltonian_] := Eigensystem[hamiltonian][[2]];

(* Ground State Construction *)
negenergy[hamiltonian_] := Select[energy[hamiltonian], Negative];
negvects[hamiltonian_] := Pick[vects[hamiltonian], Negative[energy[hamiltonian]]];
groundstate[hamiltonian_, sites_] := Sqrt[2] Total[negvects[hamiltonian]]/Sqrt[sites];

(* Position Operator *)
carpos[sites_, hamiltonian_] :=  Inverse[vects[hamiltonian]] DiagonalMatrix[aa*Range[0, (sites - 1), 1]] vects[hamiltonian];

(* Polarization Operator *)
polarization[sites_,hamiltonian_] := (charge/(2 sites*aa))* ( Sum[i*aa, {i, 0, sites - 1}] - (Conjugate[groundstate[hamiltonian,  sites]].  (DiagonalMatrix[aa*Range[0, (sites - 1), 1]] .groundstate[hamiltonian,sites]) )) ;

(* Gaussian *)
gauss[\[Alpha]_,\[Epsilon]_,e_] := ( \[Alpha] / \[Pi] ) Exp[ - \[Alpha] (e - \[Epsilon])^2]; (* \[Alpha] inverse of the lenght, \[Epsilon] center of the gaussian, e variable *)

(* Density of State *)
dens[\[Alpha]_, \[Epsilon]_, energy_, sites_] :=  2/(sites*aa) Total[ gauss[\[Alpha], \[Epsilon], energy]];



(* PLOTS *)

(* Polarization per Unit Cell *)
ListPlot[Table[polarization[i, hamchainOBCdisorder[i]/(0.5*i), {i, 2, 100,2}], AxesOrigin -> {0, 0}, Frame -> True, FrameLabel -> {Style["Number of Unit Cells"], Style["Polarization per Unit Cell"]}]

(* Density of States *)
Manipulate[{dos = Total[gauss[\[Alpha], energy[hamchainOBCtManipulate[{dos = Total[gauss[\[Alpha], energy[hamchainOBCdisorder[sites]], x]]; Plot[dos, {x, -5, 5}, ImageSize -> {400, 400}]}, {{\[Alpha], 2, "\[Alpha]"}, 1, 50, 1, Appearance -> "Labeled"}, {{sites, 2, "sites"}, 2, 200, 2, Appearance -> "Labeled"}, TrackedSymbols :> {\[Alpha], sites}]dep[sites, t]], x]]/sites; Plot[dos, {x, -5, 5}, ImageSize -> {400, 400}]}, {{sites, 2, "N"}, 2, 100, 2, Appearance -> "Labeled"}, {{\[Alpha], 2, "\[Alpha]"}, 1,10, 0.5, Appearance -> "Labeled"}, {{t, 0, "t"}, 0, -1, -0.01,Appearance -> "Labeled"}, TrackedSymbols :> {sites, \[Alpha], t}]
