- [ ] Usare codifica in binario.
- [ ] Fare una funzione di check che verifica che ogni riga verifica le condizioni poste.
- [ ] Una volta che si verifica il check spegnere il check
- [ ] Inventarsi operatore di selezione, scegliendo una metrica (distanza più piccola o più grande).
      L'indice scelto j è la parte intera di M * rand(0,1) ^ p.
      - Se si usa una potenza maggiore di uno, ho una legge che sale con concavità verso l'alto, a maggior parte delle selezioni fatte ha valore verso il basso, buona scelta se migliore è in fondo alla sequenza,
      - Se si usa una potenza minore di uno, spesso otterrò quelli in cima alla lista, quindi in questo caso è meglio mettere il migliore in cima alla lista.
- In base alla codifica scelta si ha una degenerazione, e quasi tutte ce l'hanno. Il miglior percorso per esempio può essere scelto in senso orario o antiorario.
    Inoltre il miglior percorso può essere iniziato da una qualunque città, quindi __DEVO FISSARE LA CITTà di PARTENZA__. Se no si rischia una soluzione antagonista ed è pericolosissiam.
- [ ] Operazioni di mutazione, vista la degenerazione e la difficoltà del problema combinatorio usa:
    - Operatore di permutazione di due città. __NON DEVE LAVORARE SULLA PRIMA CITTà.__ Posso applicarlo ripetutamente sulla prima sequenza ordinata per creare la popolazione
    - Operatore di shift di città contigue. Shifta di n posizione m città contigue. Questo può preservare una sequenza buona a differenza del precedente che è distruttivo.
    - Operature di permutazione di n città contigue. Questo può risolvere la degenerazione del senso, risolvendo mezzo percorso in un senso e mezzo nell'altro.
    - Operatore di inversione di n città contigue.
    Tieni probabilità di mutazione basse.
