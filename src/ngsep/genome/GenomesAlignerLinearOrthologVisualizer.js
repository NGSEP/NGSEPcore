const option = d3.selectAll('#option');
const form = option.append('form');
const containerSection = option.append('div').attr('class', 'container section');
containerSection.append('h1').attr('class', 'black-text center').text('Linear visualization');
option.append('input').attr('type', 'button').attr('onClick', 'history.go(0)').attr('value', 'Start again!');
option.append('input').attr('id', 'SyntenyButton').attr('type', 'button').attr('value', 'Synteny orthologs');
d3.select('#SyntenyButton').on('click', showSynteny);
option.append('input').attr('id', 'AllButton').attr('type', 'button').attr('value', 'All orthologs');
d3.select('#AllButton').on('click', showAll);

// TODO?: Slider instead of text input for better user feedback
option.append('input').attr('id', 'MinimumChromosomeLengthTextInput').attr('type', 'text').attr('placeholder', `100000`);
option.append('input').attr('id', 'MinimumChromosomeLengthButton').attr('type', 'button').attr('value', `ChromosomeLength`);
d3.select('#MinimumChromosomeLengthButton').on('click', chromosomeLength);


const canvasContainer = option.append('div').attr('class', 'container');
canvasContainer.append('div').attr('class', 'canvas');


minimumChromosomeLength = 100000;
const dims = {
    height: 800,
    width: 1000
}

const margin = { left: 200, right: 200, top: 20, bottom: 20 };

let allOrthologs = {};

const graph = d3.selectAll('.canvas')
    .append('svg')
    .attr('width', dims.width + margin.left + margin.right)
    .attr('height', dims.height + margin.top + margin.bottom)
    .attr('transform', `translate(${margin.left}, ${margin.top})`);

// group for chromosome labels
const chromosomeLabelsG1 = graph.append('g')
    .attr('class', 'chromosomeLabelsG1');
const chromosomeLabelsG2 = graph.append('g')
    .attr('class', 'chromosomeLabelsG2');

// group for synteny blocks
const blocksGroup = graph.append('g')
    .attr('class', 'syntenyBlocks');

// group for ortholog lines
const linesGroup = graph.append('g')
    .attr('class', 'orthologLines');



//scales
const y1 = d3.scaleLinear()
    .range([0, dims.height]);
const y2 = d3.scaleLinear()
    .range([0, dims.height]);
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

const kiloFormat = d3.formatPrefix(',.1', 1e3);
const megaFormat = d3.formatPrefix(',.1', 1e6);
const gigaFormat = d3.formatPrefix(',.1', 1e9);

// axes groups
const y1AxisGroup = graph.append('g')
    .attr('class', 'y1-axis')
    .attr('transform', `translate(${margin.left}, 0)`);
const y2AxisGroup = graph.append('g')
    .attr('class', 'y2-axis')
    .attr('transform', `translate(${dims.width}, 0)`);

// brush (selection)
const brush1 = d3.brushY()
    .extent([[0, 0], [dims.width / 2, dims.height]])
    .on('end', brushended1);
const brush2 = d3.brushY()
    .extent([[dims.width / 2, 0], [dims.width, dims.height]])
    .on('end', brushended2);
var idleTimeout;
const idleDelay = 350;

graph.append('g')
    .attr('class', 'brush1')
    .attr('transform', `translate(${margin.left}, 0)`)
    .call(brush1);

graph.append('g')
    .attr('class', 'brush2')
    .call(brush2);

function brushended1() {
    var s = d3.event.selection;
    if (!s) {
        if (!idleTimeout) return idleTimeout = setTimeout(idled, idleDelay);
        y1.domain([0, maxG1]);
    } else {
        y1.domain([s[0], s[1]].map(y1.invert, y1));
        graph.select('.brush1').call(brush1.move, null);
    }
    zoom1();
}

function brushended2() {
    var s = d3.event.selection;
    if (!s) {
        if (!idleTimeout) return idleTimeout = setTimeout(idled, idleDelay);
        y2.domain([0, maxG2]);
    } else {
        y2.domain([s[0], s[1]].map(y2.invert, y2));
        graph.select('.brush2').call(brush2.move, null);
    }
    zoom2();
}

function idled() {
    idleTimeout = null;
}

function zoom1() {
    var t = graph.transition().duration(750);
    if (y1.domain()[0] === 0 && y1.domain()[1] === maxG1) {
        restartTicksG1(t, false);
    } else {
        restartTicksG1(t, true);
        newTicks = [];
        let i;
        let newTickValue;
        let formatConstraint;
        for (i = 1; i <= 10; i++) {
            newTickValue = Math.round(y1.domain()[0] + (y1.domain()[1] - y1.domain()[0]) / 10 * i);
            formatConstraint = y1.domain()[1] - y1.domain()[0];
            if (newTicks.indexOf(newTickValue) === -1) {
                newTicks.push(newTickValue);
            }
        }
        y1Axis.tickValues([...y1Axis.tickValues(), ...newTicks]).tickFormat((d, i) => {
            realTickValue = d - lengthsG1[topPinnedLabelG1(d)];
            return formatConstraint > 1000000000 ? gigaFormat(realTickValue) : formatConstraint > 1000000 ? megaFormat(realTickValue) : formatConstraint > 1000 ? kiloFormat(realTickValue) : d3.format(',')(realTickValue);
        });
    }
    y1AxisGroup.transition(t).call(y1Axis);

	paintSyntenyBlocks (syntenyBlocks)
    graph.selectAll('line.orthologLine').transition(t)
        .attr('y1', d => {
            return y1(d.geneStart)
        })
        .attr('opacity', d => d.geneStart > y1.domain()[1] || d.geneStart < y1.domain()[0] || d.geneStartG2 > y2.domain()[1] || d.geneStartG2 < y2.domain()[0] ? 0.0 : 1.0);
    // .attr('visibility', d => d.geneStart > y1.domain()[1] || d.geneStart < y1.domain()[0] || d.geneStartG2 > y2.domain()[1] || d.geneStartG2 < y2.domain()[0] ? 'hidden' : null);
    topPinnedLabel = topPinnedLabelG1(y1.domain()[0]);
    bottomPinnedLabel = bottomPinnedLabelG1(y1.domain()[1]);
    graph.selectAll('text.chromosomeLabelG1').transition(t)
        .attr('transform', d => `translate(${margin.left - 40}, 
        ${d.Name === topPinnedLabel ?
                y1.range()[0] + margin.top :
                d.Name === bottomPinnedLabel ?
                    y1(lengthsG1[d.Name] + 1) + 10 :
                    y1(d.Length / 2 + lengthsG1[d.Name]) + 10})`
        );
}

function zoom2() {
    var t = graph.transition().duration(750);
    if (y2.domain()[0] === 0 && y2.domain()[1] === maxG2) {
        restartTicksG2(t, false);
    } else {
        restartTicksG2(t, true);
        newTicks = [];
        let i;
        let newTickValue;
        let formatConstraint;
        for (i = 1; i <= 10; i++) {
            newTickValue = Math.round(y2.domain()[0] + (y2.domain()[1] - y2.domain()[0]) / 10 * i);
            formatConstraint = y2.domain()[1] - y2.domain()[0];
            if (newTicks.indexOf(newTickValue) === -1) {
                newTicks.push(newTickValue);
            }
        }
        y2Axis.tickValues([...y2Axis.tickValues(), ...newTicks]).tickFormat((d, i) => {
            realTickValue = d - lengthsG2[topPinnedLabelG2(d)];
            return formatConstraint > 1000000000 ? gigaFormat(realTickValue) : formatConstraint > 1000000 ? megaFormat(realTickValue) : formatConstraint > 1000 ? kiloFormat(realTickValue) : d3.format(',')(realTickValue);
        });
    }
    y2AxisGroup.transition(t).call(y2Axis);

	paintSyntenyBlocks (syntenyBlocks)
    graph.selectAll('line.orthologLine').transition(t)
        .attr('y2', d => {
            return y2(d.geneStartG2)
        })
        .attr('opacity', d => d.geneStart > y1.domain()[1] || d.geneStart < y1.domain()[0] || d.geneStartG2 > y2.domain()[1] || d.geneStartG2 < y2.domain()[0] ? 0.0 : 1.0);
    // .attr('visibility', d => d.geneStart > y1.domain()[1] || d.geneStart < y1.domain()[0] || d.geneStartG2 > y2.domain()[1] || d.geneStartG2 < y2.domain()[0] ? 'hidden' : null);
    topPinnedLabel = topPinnedLabelG2(y2.domain()[0]);
    bottomPinnedLabel = bottomPinnedLabelG2(y2.domain()[1]);
    graph.selectAll('text.chromosomeLabelG2').transition(t)
        .attr('transform', d => `translate(${dims.width + 40}, 
            ${d.Name === topPinnedLabel ?
                y2.range()[0] + margin.top :
                d.Name === bottomPinnedLabel ?
                    y2(lengthsG2[d.Name] + 1) + 10 :
                    y2(d.Length / 2 + lengthsG2[d.Name]) + 10})`);
}

//auxiliary function to get which label to pin
const topPinnedLabelG1 = value => {
    let pinnedLabel;
    lengthsG1Map.forEach((value1,chromosome) => {
        if (value1 < value) {
            pinnedLabel = chromosome
        }
    });
    return pinnedLabel;
}
const bottomPinnedLabelG1 = value => {
    let pinnedLabel = 0;
    Object.keys(lengthsG1).forEach((chromosome, index, array) => {
        const invertedChromosome = array[array.length - index - 1];
        if (lengthsG1[invertedChromosome] < value && pinnedLabel === 0) {
            pinnedLabel = invertedChromosome
        }
    });
    if (y1.domain()[0] === 0 && y1.domain()[1] === maxG1) {
        return;
    }
    return pinnedLabel;
}

const topPinnedLabelG2 = value => {
    let pinnedLabel;
    lengthsG2Map.forEach((value2,chromosome) => {
        if (value2 < value) {
            pinnedLabel = chromosome
        }
    });
    return pinnedLabel;
}
const bottomPinnedLabelG2 = value => {
    let pinnedLabel = 0;
    Object.keys(lengthsG2).forEach((chromosome, index, array) => {
        const invertedChromosome = array[array.length - index - 1];
        if (lengthsG2[invertedChromosome] < value && pinnedLabel === 0) {
            pinnedLabel = invertedChromosome
        }
    });
    if (y2.domain()[0] === 0 && y2.domain()[1] === maxG2) {
        return;
    }
    return pinnedLabel;
}

// Create axes
// TODO: ticks need to restart at every chromosome, 
// or appear only at a certain zoom level
const y1Axis = d3.axisLeft(y1)
    .ticks(4);
const y2Axis = d3.axisRight(y2)
    .ticks(4);

let maxG1;
let maxG2;
let lengthsG1;
let lengthsG1Map;
let lengthsG2;
let lengthsG2Map;
let ticksG1;
let ticksG2;
const prepareData = () => {
    // Filter chromosomes/scaffolds
    genomeData1 = genome1.filter(chromosome => {
        return chromosome.Length > minimumChromosomeLength;
    });
    genomeData2 = genome2.filter(chromosome => {
        return chromosome.Length > minimumChromosomeLength;
    });
    
    // Get the max and relative lengths for axes and ticks
    maxG1 = 0;
    lengthsG1 = {};
	lengthsG1Map = new Map();
    ticksG1 = [];
    genomeData1.forEach(g => {
        ticksG1.push(maxG1 + 1);
        lengthsG1[g.Name] = maxG1;
		lengthsG1Map.set(g.Name, maxG1);
        maxG1 += parseInt(g.Length);
        ticksG1.push(maxG1);
    });
    
    maxG2 = 0;
    lengthsG2 = {};
	lengthsG2Map = new Map();
    ticksG2 = [];
    genomeData2.forEach(g => {
        ticksG2.push(maxG2 + 1);
        lengthsG2[g.Name] = maxG2;
		lengthsG2Map.set(g.Name, maxG2);
        maxG2 += parseInt(g.Length);
        ticksG2.push(maxG2);
    });
    // Update scale domains with the newfound max
    y1.domain([0, maxG1]);
    y2.domain([0, maxG2]);
    // Update color domain with chromosome names
    color.domain(genomeData1.map(chromosome => chromosome.Name));

    // Call axes, create ticks
    restartTicksG1(null, false);
    restartTicksG2(null, false);

    // Chromosome labels
    const labelsG1 = chromosomeLabelsG1.selectAll('text').data(genomeData1);
    labelsG1.exit().remove();
    labelsG1.enter()
        .append('text')
        .attr('class', 'chromosomeLabelG1')
        .text(d => d.Name)
        .attr('transform', d => `translate(${margin.left - 50}, ${y1(d.Length / 2 + lengthsG1[d.Name]) + 10})`)
        .attr('text-anchor', 'end')
        .attr('fill', d => color(d.Name));
    const labelsG2 = chromosomeLabelsG2.selectAll('text').data(genomeData2);
    labelsG2.exit().remove();
    labelsG2.enter()
        .append('text')
        .attr('class', 'chromosomeLabelG2')
        .text(d => d.Name)
        .attr('transform', d => `translate(${dims.width + 50}, ${y2(d.Length / 2 + lengthsG2[d.Name]) + 10})`)
        .attr('text-anchor', 'start')
        .attr('fill', d => color(d.Name));
}

const restartTicksG1 = (t, zoomed = false) => {
    y1Axis.tickValues(ticksG1).tickFormat((d, i) => {
        realTickValue = d - lengthsG1[topPinnedLabelG1(d)];
        return realTickValue > 1000000000 ? gigaFormat(realTickValue) : realTickValue > 1000000 ? megaFormat(realTickValue) : realTickValue > 1000 ? kiloFormat(realTickValue) : realTickValue != 1 ? realTickValue : zoomed ? realTickValue : null;
    });
    y1AxisGroup.transition(t).call(y1Axis);
}

const restartTicksG2 = (t, zoomed = false) => {
    y2Axis.tickValues(ticksG2).tickFormat((d, i) => {
        realTickValue = d - lengthsG2[topPinnedLabelG2(d)];
        return realTickValue > 1000000000 ? gigaFormat(realTickValue) : realTickValue > 1000000 ? megaFormat(realTickValue) : realTickValue > 1000 ? kiloFormat(realTickValue) : realTickValue != 1 ? realTickValue : zoomed ? realTickValue : null;
    });
    y2AxisGroup.transition(t).call(y2Axis);
}

var lineFunc = d3.line()
  .x(function(d) { return d.x })
  .y(function(d) { return d.y })

const paintSyntenyBlocks = blocks => {
	// Create data to draw the lines in the correct positions
	blocksData = []
    blocks.forEach(block => {
    	    
        if (chromosomesDisplayed(genomeData1, genomeData2, block.chromosomeG1, block.chromosomeG2)) {
			y11 = lengthsG1[block.chromosomeG1] + parseInt(block.regionStartG1);
			y12 = lengthsG1[block.chromosomeG1] + parseInt(block.regionEndG1);
			y21 = lengthsG2[block.chromosomeG2] + parseInt(block.regionStartG2);
			y22 = lengthsG2[block.chromosomeG2] + parseInt(block.regionEndG2);
			if(y11>=y1.domain()[0] && y12 <=y1.domain()[1] && y21>=y2.domain()[0] && y22 <=y2.domain()[1]) {
				block.coords = [];
                block.coords.push( { 'x': margin.left, 'y': y1(y11) } );
				if (block.negativeG2==1) {
					block.coords.push( { 'x': dims.width, 'y': y2(y22) } );
					block.coords.push( { 'x': dims.width, 'y': y2(y21) } );
				} else {
					block.coords.push( { 'x': dims.width, 'y': y2(y21) } );
                	block.coords.push( { 'x': dims.width, 'y': y2(y22) } );
				}
                block.coords.push( { 'x': margin.left, 'y': y1(y12) } );
                blocksData.push(block);
			}
	        
        }
    });
	// Remove previous blocks
	blocksGroup.selectAll('path.syntenyBlock').remove();
	
    // Create data to be painted
    const blockPaths = blocksGroup.selectAll('path.syntenyBlock').data(blocksData);

    // Exit selection
    blockPaths.exit().remove();

    blockPaths.enter()
        .append('path')
        .attr('d', d => lineFunc(d.coords))
        .attr('class', 'syntenyBlock')
        .attr('stroke', d => color(d.chromosomeG1))
        .attr('fill', d => color(d.chromosomeG1));
};

const paintData = orthologs => {
    // Create data to be painted
    const lines = linesGroup.selectAll('line.orthologLine').data(orthologs);

    // Exit selection
    lines.exit().remove();

    var t = graph.transition().duration(750);

    lines
        // .transition(t) // Uncomment this line for cool and inefficient animations
        .attr('x1', margin.left)
        .attr('x2', dims.width)
        .attr('y1', d => y1(d.geneStart))
        .attr('y2', d => y2(d.geneStartG2))
        .style('stroke', d => color(d.chromosome));

    lines.enter()
        .append('line')
        .attr('class', 'orthologLine')
        .attr('id', d => `${d.geneId}::${d.geneIdG2}`)
        .attr('x1', margin.left)
        .attr('x2', dims.width)
        .attr('y1', d => y1(d.geneStart))
        .attr('y2', d => y2(d.geneStartG2))
        .style('stroke', d => color(d.chromosome));

    graph.selectAll('line.orthologLine')
        .attr('y1', d => {
            return y1(d.geneStart)
        })
        .attr('opacity', d => d.geneStart > y1.domain()[1] || d.geneStart < y1.domain()[0] || d.geneStartG2 > y2.domain()[1] || d.geneStartG2 < y2.domain()[0] ? 0.0 : 1.0);
};

// Auxiliary function to check if a chromosome is displayed
const chromosomesDisplayed = (genomeData1, genomeData2, chrG1, chrG2) => {
    const displayed = { genome1: false, genome2: false };
    genomeData1.forEach(chromosome => {
        if (chromosome.Name == chrG1) {
            displayed.genome1 = true;
        }
    })
    genomeData2.forEach(chromosome => {
        if (chromosome.Name == chrG2) {
            displayed.genome2 = true;
        }
    })
    return displayed.genome1 && displayed.genome2;
};


// This is actually a filter for the orthologs, but can be extended
const createLineData = orthologs => {
    lineData = []
    orthologs.forEach(ortholog => {
        if (chromosomesDisplayed(genomeData1, genomeData2, ortholog.chromosome, ortholog.chromosomeG2)) {
            ortholog.geneStart = lengthsG1[ortholog.chromosome] + parseInt(ortholog.geneStart);
            ortholog.geneStartG2 = lengthsG2[ortholog.chromosomeG2] + parseInt(ortholog.geneStartG2);
            ortholog.geneEnd = lengthsG1[ortholog.chromosome] + parseInt(ortholog.geneEnd);
            ortholog.geneEndG2 = lengthsG2[ortholog.chromosomeG2] + parseInt(ortholog.geneEndG2);

            lineData.push(ortholog);
        }
    });
    return lineData;
};

// Show synteny
function showSynteny() {
	blocksGroup.selectAll('path.syntenyBlock').remove();
	paintData(allOrthologs.synteny);
	
    document.getElementById('SyntenyButton').disabled = true;
    document.getElementById('AllButton').disabled = false;
}

// Show all
function showAll() {
	blocksGroup.selectAll('path.syntenyBlock').remove();
    paintData(allOrthologs.all);
    document.getElementById('SyntenyButton').disabled = false;
    document.getElementById('AllButton').disabled = true;
}

// Filter chromosomes
function chromosomeLength() {
    const newMinimum = d3.select('#MinimumChromosomeLengthTextInput').property('value');
    if (newMinimum != '') {
        minimumChromosomeLength = parseInt(newMinimum);
        prepareData();
    }
}

// Split data
const divideOrthologs = (orthologs,genomeId2) => {
    orthologsSynteny = [];
    orthologsAll = [];
    orthologs.map(ortholog => {
	if (parseInt(ortholog.genomeId2) == genomeId2) {
	    if (ortholog.block > -1) {
	        orthologsSynteny.push(ortholog)
	    }
	    orthologsAll.push(ortholog);
	}
    });
    orthologs = createLineData(orthologsAll);

    dividedOrthologs = { 'synteny': orthologsSynteny, 'all': orthologsAll };
    return dividedOrthologs;
}

// Read data
prepareData();
allOrthologs = divideOrthologs(orthologsG1,2);
paintSyntenyBlocks(syntenyBlocks);
