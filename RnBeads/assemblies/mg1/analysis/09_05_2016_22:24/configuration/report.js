function toggleSection(sid) {
	var elh = document.getElementById('hsection' + sid)
	var eld = document.getElementById('section' + sid)
	if (elh.className == "collapsed") {
		elh.className = "expanded"
		eld.style.display = "block"
	} else { // el.className == "expanded"
		elh.className = "collapsed"
		eld.style.display = "none"
	}
}

function doCursor() {
	if (document.body.style.cursor == 'pointer') {
		document.body.style.cursor = 'auto'
	} else {
		document.body.style.cursor = 'pointer'
	}
}

function showHide(elementid) {
	var el = document.getElementById(elementid)
	if (el.style.display == 'block') {
		el.style.display = 'none'
	} else {
		el.style.display = 'block'
	}
}

function pdficon(el, activate) {
	var i = el.src.lastIndexOf('/')
	var imgprefix = (i == -1) ? '' : (el.src.substr(0, i) + '/')
	var imgsuffix = (typeof activate === 'undefined' || activate) ? 'active' : 'inactive'
	el.src = imgprefix + 'pdf' + '_' + imgsuffix + '.png'
}

function updateFigure(figid) {
	var imageel = document.getElementById(figid + 'image')
	var imagesrc = imageel.src
	var i = imagesrc.lastIndexOf('/')
	var imaged = ''
	var imagef = ''
	if (i != -1) {
		imaged = imagesrc.substr(0, i)
		imagef = imagesrc.substr(i + 1)
	}
	i = imagef.lastIndexOf('.')
	var imageext = imagef.substr(i)
	imagef = imagef.substr(0, i).split('_')
	for (i = 0; i < imagef.length; i++) {
		var selectel = document.getElementById(figid + 'setting' + (i + 1))
		if (selectel) {
			imagef[i] = selectel.value
		}
	}
	imagef = imagef.join('_')
	imageel.src = imaged + '/' + imagef + imageext
	var imagep = document.getElementById(figid + 'pdf')
	if (imagep) {
		i = imagep.href.lastIndexOf('/')
		imaged = (i == -1) ? '' : imagep.href.substr(0, i + 1)
		imagep.href = imaged + imagef + '.pdf'
	}
	var imageh = document.getElementById(figid + 'imagehigh')
	if (imageh) {
		i = imageh.href.lastIndexOf('/')
		imaged = (i == -1) ? '' : imageh.href.substr(0, i + 1)
		imagef = imaged + imagef
		if (imagef + imageext == imageel.src) {
			imagef = imagef + '_high_resolution'
		}
		imageh.href = imagef + imageext
	}
}

function updateTable(tabid) {
	var tabcandidates = document.getElementById(tabid + 'figure').getElementsByTagName('div')
	if (tabcandidates.length == 0) return

	var tableel = tabcandidates[0].getAttribute('id').substr(tabid.length + 1).split('_')
	for (i = 0; i < tableel.length; i++) {
		var selectel = document.getElementById(tabid + 'setting' + (i + 1))
		if (selectel) {
			tableel[i] = selectel.value
		}
	}

	var visibleid = tabid + '_' + tableel.join('_')
	for (i = 0; i < tabcandidates.length; i++) {
		if (tabcandidates[i].getAttribute('id') == visibleid) {
			tabcandidates[i].style.display = 'block'
		} else if (tabcandidates[i].style.display != 'none' && tabcandidates[i].style.display != '') {
			tabcandidates[i].style.display = 'none'
		}
	}
}
