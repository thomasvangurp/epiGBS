import urllib,urllib2

url = 'http://www.uniprot.org/uniparc/'

params = {
'from':'P_GI',
'to':'ACC',
'format':'tab',
'query':'590650093'
}

data = urllib.urlencode(params)
request = urllib2.Request(url, data)
contact = "" # Please set your email address here to help us debug in case of problems.
request.add_header('User-Agent', 'Python %s' % contact)
response = urllib2.urlopen(request)
page = response.read(200000)
print ''