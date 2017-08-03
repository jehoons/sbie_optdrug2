### Table S1. Fumia network analysis 

#### (A) Fumia network preparing

*A1. fumia model maually curated*

*A2. fumia model processed weight sum*

*A3. fumia model processed logical*

#### (B) Attractor analysis

*B1. Simulation*

*B2. Simulation summary*

페노타입을 결정하는 규칙은 다음과 같다: 

```
if attr_type == 'cyclic': 
    if correct_sequence: 
        P
    else: 
        Q
else:
    if apoptosis == True: 
        A
    else: 
        Unknown 
```

> 발견된 문제 - 논문에 의하면 00000은 Q상태로 100% 가지만, 시뮬레이션에서는 P상태로 100%가 되는 문제가 있다. See Figure 1. 


#### (B) Attractor analysis - *Update 1*

*B1. Simulation*

*B2. Simulation summary*

> *Ver 1*의 문제는 거의 해결된 것으로 보인다. 차이점은, 아래 입력 템플릿에서 `State_Gli: False` 로 세팅한 것이다. 그리고 point attractor에서 apoptosis 값이 false인 경우, Cyclin A,B,D,E값의 상태가 모두 False라면 Q페노타입을 가진다고 설정하였다.

```
Input = {
    'State_Mutagen' : False,
    'State_GFs': False,
    'State_Nutrients': False,
    'State_TNFalpha': False,
    'State_Hypoxia': False,
    'State_Gli': False,
    }
```

페노타입을 결정하는 규칙은 다음과 같다: 

```
if attr_type == 'cyclic': 
    if correct_sequence: 
        P
    else: 
        Q
else:
    if apoptosis == True: 
        A
    else: 
        if cyc_abde == False: 
            Q
        else: 
            Unknown 

```

#### (C) Mutational sequence analysis 

*C1. Simulation with Reduced number of combinations*

Use `get_mutations_config()`

*C2. Summary for C1*

*C3. All possible combinations*

Use `get_mutations_config_all()` 

*C4. Summary for C3*






