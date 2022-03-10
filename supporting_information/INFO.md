## Supplemetary Figures

`S1_Fig.pdf` **Random walk path length.**

`S2_Fig.pdf` **Scree plot.**

`S3_Fig.pdf` **Neurotransmitter pathways.**

`S4_Fig.pdf` **3D spatial embedding of each class.**
<!--- Each connectivity-class is visualed in a 3D template of the Drosophila brain using the natverse package. For each class we show the embeddings of each of its constieunt neruons. For each class we registered the Each neuron. The dominant neuropils for each class are also highlighted (highlighted (as a yellow colored underlay)  The Virtual Fly Brain interface allows users to explore the structure of the by browsing 3D images of a brain with subregions displayed as coloured overlays.-->

## Excel Files

`S1_File.xlsx` **ARI results**. Pairwise ARI values for different combinations of parameter choices.

`S2_File.xlsx` : **Class labels and data**
  - Sheet *Labels*: Assigned descriptor labels for each connectivity-based class.
  - Sheet *Label Description*: Brief description of the methods used to determine these labels
  - Sheet *data_conditional_prob*: The posterior probability $p(y_j|x_i)$, e.g., that the neuron is associated with the j-th functional-anatomical community, given that it belongs to the i-th connectivity-based class. Also, lists the prior probabilites for comparison.
  - Sheet *data_graph_theroetic_measures*: Multiple graph theoretic measures calculated for each connecitiy-based class. Also, lists the pairwise Absorption and Driftiness values claculated using random walks on the mesoscopic circuit.
  - Sheet *data_spad*: Statistics for the neuropil (spatial) distribution of each class.
  - Sheet *data_morph+pvec*: Statistics for the morphological and persistence vector of each class.
  - The red cells highlight the column minimum. While, the green cells highlight the column maximum.
<!--- Additional sheets include supplemental data calculated for each class, along with a summary of the procedure used to assign the corresponding labels. -->


## Movies

`S1_Movie.mp4` **Critical growth per day.**: The clip consists of 10 circuits, each corresponding to one day. This clip tracks the growth (and activity) of the circuit over the span of 10 days. If a node in the circuit has no birthed neurons until that particular day, then that node is "grayed out" and considered to be as yet inactive (unborn). For each day, the nodes that exhibit a critical growth period are colored, and so are the critical edges.

`S2_Movie.mp4` **Critical growth time periods.**: Similar to the S1_Movie but instead of each day, it shows circuit growth grouped into four time periods, (i) embryo, (ii) early (day1-day3), (iii) mid (day4-day6), and (iv) late (day7-day9).


