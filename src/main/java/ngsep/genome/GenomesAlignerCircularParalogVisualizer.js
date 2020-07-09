const option = d3.selectAll("#option");
const containerSection = option.append('div').attr('class', 'container section');
containerSection.append('h5').attr('class', 'black-text center').text('Chord visualization');
const canvasContainer = option.append('div').attr('class', 'container');
canvasContainer.append('div').attr('class', 'canvas');

const dims = {
    height: 900,
    width: 900,
    innerRadius: 270,
    outerRadius: 290,
    labelRadius: 270 + 120,
    padAngle: 0.02,
    ticksPadAngle: 0.014,
    ticksSpacing: 200000,
    ribbonPadAngle: 0.015,
    opacity: 0.7,
    fadedOpacity: 0.01,
    focusOpacity: 0.9
};

const minimumChromosomeLength = 500000;

const cent = {
    x: (dims.width / 2 + 5),
    y: (dims.height / 2 + 5)
}

const svg = d3.selectAll('.canvas')
    .append('svg')
    .attr('width', dims.width)
    .attr('height', dims.height);

const graph = svg.append('g')
    .attr('transform', `translate(${cent.x}, ${cent.y})`);

const circularAxisGroup = graph.append('g').attr('class', 'axis');

const arcsGroup = graph.append('g').attr('class', 'chromosomeArcs');

const ribbonsGroup = graph.append('g').attr('class', 'ribbons');

const pie = d3.pie()
    .sort(null)
    .value(d => d.Length);

const arc = d3.arc()
    .innerRadius(dims.innerRadius)
    .outerRadius(dims.outerRadius)
    .padAngle(dims.padAngle);

const ribbon = d3.ribbon()
    .radius(dims.innerRadius);

const color = d3.scaleOrdinal([
    '#66c2a5',
    '#fc8d62',
    '#8da0cb',
    '#e78ac3',
    '#8dd3c7',
    '#bebada',
    '#fb8072',
    '#80b1d3',
    '#fdb462',
    '#b3de69',
    '#fccde5',
    '#d9d9d9'
]);

const kiloFormat = d3.formatPrefix(',.0', 1e3);
const megaFormat = d3.formatPrefix(',.0', 1e6);
const gigaFormat = d3.formatPrefix(',.0', 1e9);

const update = (genomeData, paralogsData) => {
    genomeData = genomeData.filter(chromosome => {
        return chromosome.Length > minimumChromosomeLength;
    });
    const chordGroupData = pie(genomeData);
    let chordTicks = [];
    chordGroupData.map(item => createTicks(item, dims.ticksSpacing)).forEach(item => {
        chordTicks = [...chordTicks, ...item];
    });
    
    ribbonsData = createParalogChords(genomeData, paralogsData, chordGroupData);
    
    const ribbons = ribbonsGroup.selectAll('.ribbon').data(ribbonsData);
    const arcs = arcsGroup.selectAll('.arc').data(chordGroupData);
    const axis = circularAxisGroup.selectAll('g').data(chordTicks);

    // Update scales domains
    color.domain(chordGroupData.map(chordGroupMember => chordGroupMember.index));

    // Exit selection
    ribbons.exit().remove();
    arcs.exit().remove();

    // Current selection
    ribbons.remove();
    arcs.remove();

    // Enter selection
    const arcsToDraw = arcs.enter()

    arcsToDraw
        .append('g')
        .append('path')
        .attr('fill', d => color(d.index))
        .attr('stroke', d => d3.color(color(d.index)).darker(1))
        .style('opacity', dims.opacity)
        .attr('class', 'arc')
        .attr('id', (d, i) => `arc${i}`)
        .attr('d', arc)

    const chromosomeLabels = arcsToDraw
        .append('g')
        .each(d => { d.angle = ((d.startAngle + d.endAngle) / 2) })
        .attr("class", "outer-labels")
        .attr("text-anchor", d => d.angle > Math.PI ? "end" : null)
        .attr("transform", function (d, i) {
            const c = arc.innerRadius(dims.labelRadius).centroid(d);
            return `translate(${c[0]}, ${c[1]}) rotate(${d.angle * 180 / Math.PI - 90}) ${d.angle > Math.PI ? 'rotate(180)' : ''}`;
        })
        .attr('fill', d => color(d.index))
        .style('width', 20)

    chromosomeLabels.append("text")
        .attr("class", "outer-label")
        .attr("dy", ".35em")
        .text(d => d.data.Name)

    ribbons.enter()
        .append('g')
        .append('path')
        .attr('fill', d => color(d.target.index))
        .attr('stroke', d => d3.color(color(d.target.index)).darker(1))
        .style('opacity', dims.opacity)
        .attr('class', 'ribbon')
        .attr('d', ribbon);

    ticks = axis
        .enter()
        .append('g')
        .attr('transform', d => `rotate(${(d.angle + dims.ticksPadAngle) * 180 / Math.PI - 90}) translate(${dims.outerRadius},0)`);

    ticks
        .append('line')
        .attr('class', 'tick')
        .attr('stroke', '#000')
        .attr('x2', 5);

    ticks
        .filter(d => d.value % 5000 === 0)
        .append('text')
        .attr('class', 'tickLabel')
        .attr('x', 8)
        .attr('y', 5)
        .attr('transform', d => d.angle > Math.PI ? `rotate(180) translate(-15)` : null)
        .attr('text-anchor', d => d.angle > Math.PI ? 'end' : null)
        .text(d => d.value > 1000000000 ? gigaFormat(d.value) : d.value > 1000000 ? megaFormat(d.value) : kiloFormat(d.value))


    // Animations
    ribbonsAnimation = graph.selectAll('.ribbon');
    ribbonsGroupAnimation = graph.selectAll('.ribbons');
    arcsAnimation = graph.selectAll('.arc');

    // ribbonsAnimation
    //     .on('mouseover', (d, i, n) => {
    //         setOpacity(arcsAnimation.filter(arc => arc.index === d.target.index || arc.index === d.source.index), dims.focusOpacity);
    //         setOpacity(d3.select(n[i]), dims.focusOpacity);
    //     });

    // ribbonsGroupAnimation
    //     .on('mouseover', () => {
    //         setOpacity(ribbonsAnimation, dims.fadedOpacity);
    //         setOpacity(arcsAnimation, dims.fadedOpacity);
    //     })
    //     .on('mouseleave', () => {
    //         setOpacity(ribbonsAnimation, dims.opacity);
    //         setOpacity(arcsAnimation, dims.opacity);
    //     });

    arcsAnimation
        .on('mouseover', (d, i, n) => {
            setOpacity(arcsAnimation, dims.fadedOpacity);
            setOpacity(ribbonsAnimation, dims.fadedOpacity);
            arcsToFocus = {};
            setOpacity(ribbonsAnimation
                .filter(ribbon => {
                    if (ribbon.source.index === d.index || ribbon.target.index === d.index) {
                        if (!(ribbon.source.index in arcsToFocus)) {
                            arcsToFocus[ribbon.source.index] = n[ribbon.source.index];
                        }
                        if (!(ribbon.target.index in arcsToFocus)) {
                            arcsToFocus[ribbon.target.index] = n[ribbon.target.index];
                        }
                        return true;
                    }
                    return false;
                }), dims.focusOpacity);
            setOpacity(d3.selectAll(Object.values(arcsToFocus)), dims.focusOpacity);
        })
        .on('mouseleave', () => {
            setOpacity(arcsAnimation, dims.opacity);
            setOpacity(ribbonsAnimation, dims.opacity);
        })
};


const createTicks = (d, step) => {
    const k = (d.endAngle - d.startAngle - dims.padAngle) / d.value;
    return d3.range(0, d.value, step).map(value => {
        return { value: value, angle: value * k + d.startAngle }
    });
};

// create the actual chords.
// A chord has:
//  - Start index
//  - Source chromosome starting angle
//  - Source chromosome end angle
//  - End index
//  - Target chromosome start angle
//  - Target chromosome end angle
const createParalogChords = (genomeData, paralogs, chromosomes) => {
    paralogChords = [];
    paralogs.forEach(paralog => {
        if (chromosomesDisplayed(paralog.chromosome, paralog.paralogChr, chromosomes)) {
            const chromosomeIndexes = chromosomeToIndex(genomeData);
            targetIndex = chromosomeIndexes[paralog.paralogChr];
            sourceIndex = chromosomeIndexes[paralog.chromosome];
            sourceChromosomeStartAngle = chromosomes[sourceIndex].startAngle + dims.ribbonPadAngle;
            sourceChromosomeEndAngle = chromosomes[sourceIndex].endAngle - dims.ribbonPadAngle;
            targetChromosomeStartAngle = chromosomes[targetIndex].startAngle + dims.ribbonPadAngle;
            targetChromosomeEndAngle = chromosomes[targetIndex].endAngle - dims.ribbonPadAngle;
            sourceRibbonStartAngle = paralog.geneStart / chromosomes[sourceIndex].value * (sourceChromosomeEndAngle - sourceChromosomeStartAngle) + sourceChromosomeStartAngle;
            sourceRibbonEndAngle = paralog.geneEnd / chromosomes[sourceIndex].value * (sourceChromosomeEndAngle - sourceChromosomeStartAngle) + sourceChromosomeStartAngle;
            targetRibbonStartAngle = paralog.paralogStart / chromosomes[targetIndex].value * (targetChromosomeEndAngle - targetChromosomeStartAngle) + targetChromosomeStartAngle;
            targetRibbonEndAngle = paralog.paralogEnd / chromosomes[targetIndex].value * (targetChromosomeEndAngle - targetChromosomeStartAngle) + targetChromosomeStartAngle;
            paralogChords.push({
                source: { index: sourceIndex, subIndex: targetIndex, startAngle: sourceRibbonStartAngle, endAngle: sourceRibbonEndAngle },
                target: { index: targetIndex, subIndex: sourceIndex, startAngle: targetRibbonStartAngle, endAngle: targetRibbonEndAngle }
            });
        }
    }
    );
    return paralogChords;
}

const chromosomesDisplayed = (chromosome, paralogChromosome, chromosomes) => {
    const displayed = { chromosome: false, paralogChromosome: false };
    chromosomes.forEach((chromosomeDrawn) => {
        if (chromosome === chromosomeDrawn.data.Name) {
            displayed.chromosome = true;
        }
        if (paralogChromosome === chromosomeDrawn.data.Name) {
            displayed.paralogChromosome = true;
        }
    });
    return displayed.chromosome && displayed.paralogChromosome;
};

const chromosomeToIndex = (genomeData) => {
    let chromosomeToIndex = {};
    let counter = 0;
    genomeData.forEach(chromosome => chromosomeToIndex[chromosome.Name] = counter++);
    return chromosomeToIndex;
};

const setOpacity = (elements, opacity) => {
    elements.style('opacity', opacity);
};

update(genome1, paralogsG1);
