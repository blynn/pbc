# The following git_ref should pointing at a tag, but there is none
# Meanwhile, we can use specific branch/commit from rpmbuild
%{!?_git_ref:%define _git_ref master}

Name:           libpbc
Version:        0.5.14
Release:        0.2
BuildArch:      x86_64
Summary:        Pairing-Based Cryptography Library
Group:          System Environment/Libraries
License:        ASL 2.0
Url:            https://crypto.stanford.edu/pbc/
Source0:        https://github.com/digital-me/pbc/archive/%{_git_ref}.zip

BuildRoot:      %{_tmppath}/%{name}-%{version}-build
Requires:       gmp
Provides:       libpbc0

%description
The PBC (Pairing-Based Cryptography) library is a free C library built on the GMP library that performs the mathematical operations underlying pairing-based cryptosystems.

The PBC library is designed to be the backbone of implementations of pairing-based cryptosystems, thus speed and portability are important goals. It provides routines such as elliptic curve generation, elliptic curve arithmetic and pairing computation. Thanks to the GMP library, despite being written in C, pairings times are reasonable.

%package devel
Summary:        Headers to develop code using PBC library.
Group:          Development/Libraries
Requires:       %{name} = %{version}-%{release}
Provides:       libpbc-dev
BuildArch:      noarch

%description devel
Header files for using the PBC (Pairing-Based Cryptography) library in applications.

%prep
%setup -q -n pbc-%{_git_ref}

%build
./setup
%configure
make

%install
[ "%{buildroot}" != / ] && rm -rf "%{buildroot}"
%make_install
install -dm 755 %{buildroot}

%clean
rm -rf %{buildroot}

%pre

%post
/sbin/ldconfig

%preun

%postun
/sbin/ldconfig

%files
%defattr(-,root,root,-)
%_libdir/*

%files devel
%defattr(-,root,root,-)
/usr/include/pbc

%changelog
* Tue Mar  6 2018 - benoit.donneaux (at) digital-me.nl - 0.5.14-0.2
- Based on fork from https://github.com/digital-me/pbc.git
- Support redhat packaging and fixes debian packaging

* Wed Aug 31 2016 - benoit.donneaux (at) digital-me.nl - 0.5.14-0.1
- Rollout version 0.5.14. See https://crypto.stanford.edu/pbc/news.html
