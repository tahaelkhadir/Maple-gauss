# Maple-gauss
# by cosmosse

gauss_elimination := proc (A, B)
    local n, m, WA, i, k, solution;

    # Vérifier si la matrice A est vide
    if numelems(A) = 0 then
        error "La matrice A ne peut pas être vide";
    end if;

    # Obtenir le nombre de lignes et de colonnes de la matrice A
    n := LinearAlgebra:-RowDimension(A);
    m := LinearAlgebra:-ColumnDimension(A);

    # Créer une matrice augmentée WA
    WA := Matrix(n, n+1, 0);
    WA[1 .. n, 1 .. m] := A;
    WA[1 .. n, n+1] := B;

    # Appliquer l'élimination de Gauss
    for i to n-1 do
        # Vérifier si le pivot est zéro (éviter la division par zéro)
        if WA[i, i] = 0 then
            error "Division par zéro détectée";
        end if;

        # Appliquer l'élimination
        for k from i+1 to n do
            WA[k, 1 .. n+1] := WA[k, 1 .. n+1] - WA[k, i] * WA[i, 1 .. n+1] / WA[i, i];
        end do;
    end do;

    # Vérifier si le système a une solution unique
    for i to n do
        if WA[i, i] = 0 then
            error "Le système n'a pas de solution unique";
        end if;
    end do;

    # Résoudre le système à l'aide de la substitution arrière
    solution := Vector(n, 0);
    solution[n] := WA[n, n+1] / WA[n, n];
    for i from n-1 by -1 to 1 do
        solution[i] := (WA[i, n+1] - LinearAlgebra:-DotProduct(WA[i, i+1 .. n], solution[i+1 .. n])) / WA[i, i];
    end do;

    # Retourner la solution
    return solution;
end proc;
