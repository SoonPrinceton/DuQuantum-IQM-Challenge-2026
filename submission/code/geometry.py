"""Distance-3 rotated-code geometry and check supports, extracted unchanged."""

from collections import Counter
from surface_code import make_stim_circuit

def description():
    source=make_stim_circuit(3,1,basis='Z',no_reset=True)
    coords=source.get_final_qubit_coordinates()
    data=sorted(q for q,xy in coords.items() if int(xy[0])%2 and int(xy[1])%2)
    ancillas=sorted(set(coords)-set(data))
    x_ancillas={t.value for inst in source.flattened() if inst.name=='H' for t in inst.targets_copy()}
    checks={q:set() for q in ancillas}
    for inst in source.flattened():
        if inst.name=='CX':
            qs=[t.value for t in inst.targets_copy()]
            for a,b in zip(qs[::2],qs[1::2]):
                anc,datum=(a,b) if a in checks else (b,a)
                checks[anc].add(datum)
    frequency=Counter(q for anc in x_ancillas for q in checks[anc])
    pivots={anc:min(q for q in checks[anc] if frequency[q]==1) for anc in sorted(x_ancillas)}
    cnot=[]
    for anc,pivot in pivots.items():
        cnot.append((pivot,anc))
        cnot.extend((anc,q) for q in sorted(checks[anc]-{pivot}))
        cnot.append((pivot,anc))
    assert len(cnot)==16 and len(pivots)==4
    record=[t.value for inst in source.flattened() if inst.name=='M' for t in inst.targets_copy()]
    logical=[]
    for inst in source.flattened():
        if inst.name=='OBSERVABLE_INCLUDE':logical.extend(record[len(record)+t.value] for t in inst.targets_copy())
    coords_to_anc={tuple(coords[q][:2]):q for q in ancillas}
    z_order=[coords_to_anc[tuple(xy[:2])] for xy in source.get_detector_coordinates().values() if xy[2]==0]
    assert len(z_order)==4 and all(q not in x_ancillas for q in z_order)
    return dict(data=data,ancillas=ancillas,coords=coords,x_ancillas=sorted(x_ancillas),
                checks=checks,pivots=pivots,cnot=cnot,logical_z=logical,z_order=z_order)
