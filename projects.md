# <img border="0" src="https://www.svgrepo.com/show/1025/task.svg" width="40" height="40"> Projects

***

<br/>

Since each of the 4 group themes may have completely different softwares, please follow in instructions specific to each one.

***

### Learning strategy

The proposed educational plan for this course will be done using the Project-Based Learning (PBL) approach. This way, rather than having lectures on a topic and defined set of instructions and pre-given steps to follow (which focus on passive individual learning), you are required to learn the topic in a more **active** and **dynamic** way in order to solve the data analysis tasks as a group. Therefore, each group is provided with:

1. The <img border="0" src="https://static.thenounproject.com/png/67360-200.png" width="10" height="10">[report file](project_velocity/README.md) containing a list of instructions and questions giving the overall direction and order of the steps needed to be done (read below how to use it).
2. A <img border="0" src="https://d1nhio0ox7pgb.cloudfront.net/_img/o_collection_png/green_dark_grey/512x512/plain/dictionary.png" width="10" height="10">[glossary of terms](single_cell/glossary/glossary_of_terms_single_cell.html) which contains detailed information about many of the steps and provides further references on the topics.

***

### How to work in groups

Each group is required to gather information from the **glossary** or from other **sources** to perform single cell data analysis by filling the report file (using this <img border="0" src="https://static.thenounproject.com/png/67360-200.png" width="10" height="10">[project_report.Rmd](single_cell/code/project_report.Rmd) file), by:

1. Reading the instructions/questions on the **project report** for a certain task.
2. Searching the **glossary** or **other sources** for guidance.
3. Discussing with your colleges which step needs to be done or included.
4. Dividing tasks among your group colleagues (each one can test parts of the code and then send the result to the one responsible for the report).
5. Replacing the instructions with a text explaining the rationale for this step.
6. Adding and running code from that step.

You are encouraged to use any other additional website / article in your report. That also includes addition of other previous experiences and/or additional code steps from other sources, as long as the rationale for the inclusion is also discussed in the report text. You can also add any other plots and visualisations in your report to illustrate your results.

We recommend to save the last 1 hour available to work on the project presentations. To save time, please create a google-slide presentation and divide the tasks, so everyone can contribute to building the presentation simultaneously.

***

### Tips for a good group dynamic

Working in groups is not a simple task, but it is how science is built. The main goal of working in groups is to teach and learn with others. Be kind and patient with your colleagues, help each other so the group reaches a good and fluent dynamic. However, time to work on the project is restricted, so please focus on the tasks at hand and be clear and concise as possible when communicating. Remember: this is a group work and you are evaluated as a group player.

Here we have some tips for successfully structuring your project work:

1. The group votes one person to GUIDE the group (the most experienced is preferable). The GUIDE's functions are:
  - To keep the pace of the tasks in hand.
  - Distribute the tasks to the other group members
  - Gather results from the group members to assemble the project report.
2. Then, help each other to startup `jupyter notebook`/`Rmarkdown`, then copy the project report instructions into it.
3. Once everyone is set, the GUIDE reads the first milestone instructions out loud to the group.
4. The group then discuss what needs to be done. It is worth dividing the tasks, and the GUIDE can divide.
  - Example: When asked `which dimensionality reduction to use?`, each person finds out how to run 1 different dimensionality reduction and then **send the code and explanation to the GUIDE**, so he can paste into the report.
  - Example: When asked `which metric is best to evaluate a method?`, each person chooses 1 different technique and then **send the code and explanation to the GUIDE**, so he can paste into the report.
  - In some cases, there might be many options to choose from, so need to prioritise what to focus on.
  - If something goes wrong or is not working, ask for help from your colleagues
  - You are free to discuss in the group to skip some parts for the sake time (maybe the group tested enough for that task)
  - The GUIDE can abstain from doing this part and focus on the report, if necessary.
5. If you are done before others, you can now let the group know you are done and:
  - help your colleagues by discussing their issues
  - trying to run the code from your colleague
  - run an optional task that no one is yet running
  - help the GUIDE in writing and/or polish the report text/code
  - write code to illustrate the results and/or compare results
6. Once the task is completed and attached to the report, the GUIDE can share the updated report with the group.
7. By the end of each group session, it is beneficial that group has the latest version of the report, so you can always start from the last step.

***

### Project evaluation

By the end of the course each group is expected to have:

- **1 project report file** with all analysis steps and the rationale for the following projects and
- **1 project presentation file** showing the work done and the results achieved (10 min + 5 min questions)

Those files will be made public and hosted permanently on the course webpage. You can now even look what other groups did, see their report and go through their project.

# <img border="0" src="/single-cell_sib_scilifelab_2021/logos/single_cell.png" width="40" height="40"> Spatial Transcriptomics
***

<br/>

<details>
<summary>
<b>PROJECT 1: Blocking myeloid development during colitis</b>
</summary>

<br/>

**Background:** Ulcerative colitis (UC) is an inflammatory bowel disease (IBD) driven mainly by colonic innate inflammatory cells such as macrophages, monocytes and neutrophils ( [Czarnewski et al 2019](https://www.nature.com/articles/s41467-019-10769-x), [Skatteborg et al 2020](https://academic.oup.com/ecco-jcc/advance-article-abstract/doi/10.1093/ecco-jcc/jjaa121/5859161?redirectedFrom=fulltext)). A recent study showed that patients that present higher neutrophilic inflammatory signature (known as UC1) become refractory to both anti-TNF and anti-a4b7 integrin therapy ( [Czarnewski et al 2019](https://www.nature.com/articles/s41467-019-10769-x), [Skatteborg et al 2020](https://academic.oup.com/ecco-jcc/advance-article-abstract/doi/10.1093/ecco-jcc/jjaa121/5859161?redirectedFrom=fulltext)), which leads to surgical intervention for removal of the colon. Neutrophils inflammatory cells are short lived and originate from the common myeloid progenitor (CMP) in the bone marrow and requires constant replenishment in order to sustain elevated cell number in the colon. Herein, our main goal is to identify potential gene candidates that can block either pathways of neutrophil differentiation in the bone marrow.

**Main research question:** Which genes specifically drive the differentiation of Neutrophils.

**Importance:** Identifying such genes will allow us to: 1) perform experiment in Tamoxifen-transgenic mice where those cells can be depleted during the course of colitis. 2) find potential drugs that can inhibit those genes/pathways in order to block myeloid cell differentiation during colitis in mice (with priority to already approved drugs).

**Analysis report from group 1:**
1. [project_report_colitis](single_cell/group_reports/project_report_colitis.html) ([.Rmd](single_cell/group_reports/project_report_colitis.Rmd))
2. [Presentation](single_cell/group_reports/project_presentation_colitis.pdf)

</details>

<br/>

<details>
<summary>
<b>PROJECT 2: Identifying cell and gene candidates in severe COVID-19 patients</b>
</summary>

<br/>

**Background:** COVID-19 is an infectious disease driven by the virus SARS-CoV-2, which primarily infects lung epithelial cells. However, elderly patients usually develop severe lung inflammation and lung disfunction, ultimately leading to respiratory failure ([Guan et al 2020](https://www.nejm.org/doi/full/10.1056/nejmoa2002032)). The onset of the disease is characterised by a cytokine storm comprising several inflammatory mediators ([Pedersen et al 2020](https://www.jci.org/articles/view/137647)), specially in severe cases of the disease. Many cell types orchestrate the immune response to the virus, but their relative contribution at the single-cell resolution is still unclear. Herein, our main goal is to identify which cell types and gene pathways are altered in the blood of patients with severe COVID-19.

**Main research question:** Which cell types and genes are altered when comparing blood immune cells from healthy versus COVID-19 patients.

**Importance:** Identifying such genes will allow us to: 1) better understand why severe COVID-19 patients develop stronger immune responses; 2) find potential cells for blockage or immune enhancement therapy or; 3) identify pathways that could be targeted pharmacologically.

**Analysis report from group 2:**
1. [project_report_covid19](single_cell/group_reports/project_report_covid19.html) ([.Rmd](single_cell/group_reports/project_report_covid19.Rmd))
2. [Presentation_covid19](single_cell/group_reports/project_presentation_covid19.pdf)


</details>

<br/>

# <img border="0" src="/single-cell_sib_scilifelab_2021/logos/long_read.png" width="40" height="40"> RNA velocity and trajectories
***


Please follow the [pre-course instructions](precourse.md) in order to install all necessary software and packages.

<details>
<summary>Click to expand!</summary>

  TO DO

</details>

<br/>

# <img border="0" src="/single-cell_sib_scilifelab_2021/logos/ribo_profiling.png" width="40" height="40"> Multi-omics integration
***

<details>
<summary>Click to expand!</summary>

  TO DO

</details>

<br/>

# <img border="0" src="/single-cell_sib_scilifelab_2021/logos/uv_crosslink_ip.png" width="40" height="40"> Deep learning
***

<details>
<summary>Click to expand!</summary>

  TO DO

</details>

<br/>

<br/>

### [Back to main](README.md)
