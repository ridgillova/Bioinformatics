CREATE TABLE Subject (
	SubjectID VARCHAR(15) NOT NULL, 
	Sex CHAR(1),
	Age INT,
	BMI DECIMAL(10,2),
	IR_IS CHAR(2),
	PRIMARY KEY (SubjectID)
);

CREATE TABLE Proteome (
	SubjectID VARCHAR(15) NOT NULL,
	VisitID VARCHAR(15) NOT NULL,
	Protein VARCHAR(15),
	Measurement DECIMAL(20,5),
	PRIMARY KEY (SubjectID, VisitID),
	FOREIGN KEY (SubjectID) REFERENCES Subject(SubjectID)
);

CREATE TABLE Transcriptome (
	SubjectID VARCHAR(15) NOT NULL,
	VisitID VARCHAR(15) NOT NULL,
	Transcript VARCHAR(15),
	Measurement DECIMAL(20,5),
	PRIMARY KEY (SubjectID, VisitID),
	FOREIGN KEY (SubjectID) REFERENCES Subject(SubjectID)
);

CREATE TABLE Metabolome (
	SubjectID VARCHAR(15) NOT NULL,
	VisitID VARCHAR(15) NOT NULL,
	PeakID VARCHAR(255),
	Measurement DECIMAL(20,5),
	PRIMARY KEY (SubjectID, VisitID),
	FOREIGN KEY (SubjectID) REFERENCES Subject(SubjectID)
);

CREATE TABLE Annotation (
	MetaboliteName VARCHAR(255) NOT NULL,
	KEGGID VARCHAR(15),
	Pathway VARCHAR(255),
	PRIMARY KEY (MetaboliteName)
);

CREATE TABLE Metabolite_annotation (
	PeakID VARCHAR(255) NOT NULL,
	MetaboliteName VARCHAR(255) NOT NULL,
	PRIMARY KEY (PeakID, MetaboliteName),
	FOREIGN KEY (PeakID) REFERENCES Metabolome(PeakID),
	FOREIGN KEY (MetaboliteName) REFERENCES Annotation(MetaboliteName)
);
