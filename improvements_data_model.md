# Improvements

Some thoughts derived from looking at how QUEEN does things. Trying to capture the same info with known datatypes.

* PCR: currently uses `fragment_boundaries`, `primer_footprints`. Instead, it could use two oriented features, in the same order as the primers.
* Restriction: currently uses `fragment_boundaries` and `restriction_enzymes`, instead it could use a feature for each, and add the qualifier cutsite like in this example:
    ```
    type: misc_feature
    location: [34:40](+)
    qualifiers:
        Key: cutsite, Value: ['G^TCGA_C']
    ```
* Sticky ligation: looks alright, but maybe could be made more generic, in the end maybe it may not have to be different from the description of a gibson assembly.

## Genbank features

https://www.insdc.org/submitting-standards/feature-table/

```
Feature Key           primer_bind


Definition            non-covalent primer binding site for initiation of
                      replication, transcription, or reverse transcription;
                      includes site(s) for synthetic e.g., PCR primer elements;

Optional qualifiers   /allele="text"
                      /db_xref="<database>:<identifier>"
                      /experiment="[CATEGORY:]text"
                      /gene="text"
                      /gene_synonym="text"
                      /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"
                      /locus_tag="text" (single token)
                      /map="text"
                      /note="text"
                      /old_locus_tag="text" (single token)
                      /standard_name="text"
                      /PCR_conditions="text"

Comment               used to annotate the site on a given sequence to which a primer 
                      molecule binds - not intended to represent the sequence of the
                      primer molecule itself; PCR components and reaction times may 
                      be stored under the "/PCR_conditions" qualifier; 
                      since PCR reactions most often involve pairs of primers,
                      a single primer_bind key may use the order() operator
                      with two locations, or a pair of primer_bind keys may be
                      used.
```

Note how they mention the `order()` operator. This could actually capture the orientation as well.

Have a look as well at the `source` feature key

`variation` feature key is also interesting.

## Genbank qualifiers:

* db_xref
* experiment: a brief description of the nature of the experimental evidence that supports the feature identification or assignment.
* note: comments
* PCR_conditions
* PCR_primers: <details> <summary>expand</summary>
    ```
        Qualifier       /PCR_primers=
        Definition      PCR primers that were used to amplify the sequence.
                        A single /PCR_primers qualifier should contain all the primers used  
                        for a single PCR reaction. If multiple forward or reverse primers are                   
                        present in a  single PCR reaction, multiple sets of fwd_name/fwd_seq 
                        or rev_name/rev_seq values will be  present.
        Value format    /PCR_primers="[fwd_name: XXX1, ]fwd_seq: xxxxx1,[fwd_name: XXX2,]
                        fwd_seq: xxxxx2, [rev_name: YYY1, ]rev_seq: yyyyy1, 
                        [rev_name: YYY2, ]rev_seq: yyyyy2"

        Example         /PCR_primers="fwd_name: CO1P1, fwd_seq: ttgattttttggtcayccwgaagt,
                        rev_name: CO1R4, rev_seq: ccwvytardcctarraartgttg"
                        /PCR_primers=" fwd_name: hoge1, fwd_seq: cgkgtgtatcttact, 
                        rev_name: hoge2, rev_seq: cg<i>gtgtatcttact" 
                        /PCR_primers="fwd_name: CO1P1, fwd_seq: ttgattttttggtcayccwgaagt,
                        fwd_name: CO1P2, fwd_seq: gatacacaggtcayccwgaagt, rev_name: CO1R4,  
                        rev_seq: ccwvytardcctarraartgttg" 

        Comment         fwd_seq and rev_seq are both mandatory; fwd_name and rev_name are
                        both optional. Both sequences should be presented in 5'>3' order. 
                        The sequences should be given in the IUPAC degenerate-base alphabet,
                        except for the modified bases; those must be enclosed within angle
                        brackets <> 
    ```
</details>

* standard_name: accepted standard name for this feature. Normally used to describe gene names.

