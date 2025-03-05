```mermaid
flowchart TD
    A[Start] --> B[Input: input_desk_DIC]
    B --> C[Call M_J_KIII_2D]
    C --> D1[Input: alldata, MatProp, loopedJ]
    D1 --> E{MatProp.Operation is 'DIC'?}
    E --> |Yes| F[Reshape Data]
    F --> G1[Call reshapeData(alldata)]
    G1 --> G2[Call crackgradient(RawData.Ux, stepsize)]
    G2 --> G3[Call crackgradient(RawData.Uy, stepsize)]
    G3 --> G4{alldata has 6 columns?}
    G4 --> |Yes| G5[Call crackgradient(RawData.Uz, stepsize)]
    G4 --> |No| G6[Continue]
    E --> |No| G6[Continue]
    G6 --> H[Prepare Data]
    H --> I[Check and adjust alldata format]
    I --> J[Call reshapeDefromationGradient(alldata)]
    J --> K[Set Material Properties in Maps]
    K --> L[Prompt User for Cropping and Centering Crack Tip]
    L --> M{User selects Yes?}
    M --> |Yes| N[Call CroppingEqually(Maps)]
    N --> O[Call center_Crack_tip(Crop)]
    M --> |No| P[Continue]
    L --> P[Continue]
    P --> Q[Prompt User for Crack Position]
    Q --> R{Crack on Right?}
    R --> |Yes| S[Flip Maps data]
    R --> |No| T[Continue]
    S --> T[Continue]
    T --> U[Calculate J-Integral]
    U --> V[Call decomposeDU(Maps)]
    V --> W[Calculate Work Done (Wd)]
    W --> X[Calculate Domain Integral]
    X --> Y[Contour Selection and Summation]
    Y --> Z[Calculate Stress Intensity Factors (SIFs)]
    Z --> AA[Plot and Save Results]
    AA --> AB[Call plot_JKIII(KI, KII, KIII, J, stepsize, input_unit)]
    AB --> AC[End]
```