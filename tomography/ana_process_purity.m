function purity = ana_process_purity(pm)
% function purity = ana_process_purity(pm)
% Given a process matrix pm, calculate it's purity (overlap with
% the maximum eigenvector of the Chi matrix (see http://lanl.arxiv.org/pdf/1202.5344.pdf)
c=MtoChi(pm);
[e v] = eig(c);
[~, i] = max(diag(v));
psi=e(:,i);
psi=psi/sqrt(dot(psi,psi));
% From IBM paper
purity=psi'*c*psi;
% My half-assed thing.
%purity=sqrt(sum(abs(eig(pm)).^2)/4)
end

