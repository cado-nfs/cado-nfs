import logging
import sys
import socket
import os
from subprocess import check_output, CalledProcessError, STDOUT

logger = logging.getLogger(__name__)

# Ideally, just about everything we do here should be covered by the
# python cryptography module. https://cryptography.io/ It's actually
# tempting to make it a requirement. Note that when flask itself is
# installed via pip, it does not pull the cryptography package.  It does
# not even need it to run an https server if the cert and key are
# provided separately (it does when they're requested on-the-fly with the
# ssl_context='adhoc' mode). In contrast, the debian package
# python3-flask does pull python3-cryptography)

openssl_config_template = \
r"""
oid_section		= new_oids

[ new_oids ]
[ ca ]
[ req ]
default_bits		= {bits:d}
distinguished_name	= req_distinguished_name
attributes		= req_attributes
x509_extensions	= v3_ca
string_mask = utf8only

[ req_distinguished_name ]
[ req_attributes ]
[ v3_ca ]
subjectKeyIdentifier=hash
authorityKeyIdentifier=keyid:always,issuer
basicConstraints = critical,CA:true
subjectAltName=@altnames
[altnames]
{SAN:s}
"""  # noqa: E122


def create_certificate(cafile, name, SAN, bits=2048):
    """
    calls openssl to create a certificate named cafile, with the given
    (common) name and alternate names.

    This writes the corresponding openssl config file to {cafile}.config
    and the corresponding private key to {cafile}.key. The pair of file
    names (certificate, keyfile) is returned and can be used directly as
    an ssl_context argument for Flask.run, for example (but it's really
    just a pair of bloody strings)
    """

    keyfile = cafile + '.key'

    if os.path.isfile(cafile) and os.access(cafile, os.R_OK) and \
       os.path.isfile(keyfile) and os.access(keyfile, os.R_OK):
        return (cafile, keyfile)

    config_filename = cafile + '.config'
    print(openssl_config_template.format(bits=2048, SAN=SAN),
          file=open(config_filename, 'w'))

    subject = [
        "C=XY",
        "ST=None",
        "O=None",
        "localityName=None",
        "commonName=" + name,
        "organizationalUnitName=None",
        "emailAddress=None"
    ]

    command = ['openssl', 'req', '-new', '-x509', '-batch', '-days', '365',
               '-nodes', '-subj', '/%s/' % '/'.join(subject),
               '-config', config_filename,
               '-out', cafile, '-keyout', keyfile]
    logger.debug("Running " + " ".join(command))
    try:
        output = check_output(command, stderr=STDOUT)
    except (OSError, CalledProcessError) as e:
        logger.warning(f"openssl failed: {e}")
        return None
    logger.debug("openssl output: " + output.decode('utf-8'))
    return (cafile, keyfile)


def get_server_alternate_names(address):
    """
    given a server that will listen on the given address (maybe it's an
    ipv4 address, maybe an fqdn, etc), try to determine possible server
    alt names that it would make sense to include in the certificate
    """

    # We need to find out which addresses to put as SubjectAltNames (SAN)
    # in the certificate.

    # The server address might be given by the user in one of four ways:
    # Not specified, then url_address is the (possibly short) hostname
    # Specified as short hostname
    # Specified as FQDN
    # Specified as numeric IP address

    # We can always fill in the IP address

    url_address = address if address else socket.gethostname()
    if address is None:
        address = '0.0.0.0'

    try:
        ipaddr = socket.gethostbyname(url_address)
    except socket.gaierror as e:
        logger.error("Exception trackback: %s" % e)
        logger.error("(end of Exception trackback)")
        logger.error("Cannot resolve %s -- really weird" % url_address)
        sys.exit(1)

    SAN = "IP.1 = %s\n" % ipaddr
    fqdn = socket.getfqdn(url_address)
    # If url_address was specified as IP url_address and fqdn could not find a
    # hostname for it, we don't store it. Then only the IP url_address will
    # be given in the SAN list
    dns_counter = 1
    if not fqdn == ipaddr:
        SAN += "DNS.%d = %s\n" % (dns_counter, fqdn)
        dns_counter += 1
        # If the url_address was given as a short host name, or if
        # gethostname() produced a short host name, we store that
        if url_address != fqdn and url_address != ipaddr:
            SAN += "DNS.%d = %s\n" % (dns_counter, url_address)
            dns_counter += 1
        # If localhost is given explicitly as the listen url_address, then
        # the above adds localhost to the SAN. If nothing is given, then
        # we listen on all url_addresses (including localhost), and we
        # should add "localhost" as a SAN.
        if address == "0.0.0.0":
            SAN += "DNS.%d = localhost\n" % dns_counter
            dns_counter += 1
            SAN += "IP.2 = 127.0.0.1\n"

    return SAN


def get_certificate_hash(certfile):
    if certfile is None:
        return None
    command = ['openssl', 'x509', '-in', certfile, '-fingerprint']
    try:
        output = check_output(command)
    except (OSError, CalledProcessError) as e:
        logger.error("openssl failed: %s", e)
        return None
    output_text = output.decode("ascii")
    for line in output_text.splitlines():
        if line.startswith("SHA1 Fingerprint="):
            return line.split('=', 1)[1].replace(":", "").lower()
    return None
